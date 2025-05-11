%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Response Monitoring Theta-Band Activities across Emotional Contexts in 
% Schizophrenia- and Bipolar-Spectrum Disorders
% Suzuki, Menkes, et al.
%
% Script to preprocess data
%
% Preprocessing codes/pipeline adjusted from Delorme A. EEG is better
% left alone. Sci Rep. 2023 Feb 9; 13(1):2372.
% doi: 10.1038/s41598-023-27528-0. PMID: 36759667; PMCID: PMC9911389.
%
% Completed using MATLAB 2024a & EEGLAB 2024.0, on Windows 10 Enterprise
% Fixed for minor glitch in EEGLAB eegrej.m, line 139
% https://github.com/sccn/eeglab/commit/02a306cce16716ed5db85f260271c9dde7d85893
%  Packages
%  Cleanline v2.00
%  ICLabel v1.6
%  bva-io v1.73
%  clean_rawdata v2.91
%  dipfit v5.4
%  firfilt v2.8
%
% Author: Takakuni Suzuki
% First draft on July 2023
% Last updated May 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

%%%%%%%%% Setting up the preprocessing parameters, folders, etc. %%%%%%%%%
% Paths to folders, data, and necessary programs
ParentPath = ''; % Specifcy data folder
folders = {'FARR' 'FNEG' 'FPOS'};
eeglab_path = ['']; % Specify eeglab folder
chan_locations = [eeglab_path filesep 'plugins/dipfit/standard_BESA/standard-10-5-cap385.elp']; % Channel location file

for fol = 1:length(folders)
    TaskPath = [ParentPath filesep folders{fol}];
    iPath = [TaskPath filesep 'RawData']; % Raw data location
    oPath = [TaskPath filesep 'Preprocess_Auto']; % Output location

    % Event conditions to export files
    %%%%%%%%% Event marker codes %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Stimulus markers: 1, 2
    %     1 = congruent
    %     2 = incongruent
    %
    %   Response markers: 3, 4, 9, 10
    %     3 or 4 = correct
    %     9 or 10 = error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Resp_Cor = {'S  3' 'S  4'};
    Resp_Err = {'S  9' 'S 10'};
    Resp_con = [Resp_Cor Resp_Err]; % These are the response event codes
    Stim_con = {'S  1' 'S  2'}; % These are the stimulus event codes

    epoch_window = [-1.2 1.5];
    % Notes: This is done to remove baseline and conduct a final absolute
    % threshold-based trial-by-trial rejection. A secondary aim is to remove
    % any trials that have any boundary events that we do not want included in
    % later analyses, which dictates the decision.
    % Currently set this with assumption of 3hz as lowest frequency of
    % interest: 1 sec before baseline and 1 sec after areas of interest (up
    % to 500ms).

    % EEGLAB Preprocessing Parameters
    % After several testings by TS, more liberal thresholds than Delorme
    % (2023) to balance cleaning data and retaining most trials.
    Threshold_Min_Correlation = .80; % changed from default of .85
    Threshold_ICA= .75; % changed from default of .90
    Threshold_ARS = 100; % changed from default of  20
    Threshold_Trials = 200;

    % Create a list of participants from folder names
    cd(iPath); % Read all folder names within iPath
    File_ext = '.vhdr';
    Files = dir(['*' File_ext]);
    Participants = {Files.name};

    % Create folder for non-processed EEG files minus problematic trials defined
    % as double clicks, outside of 250-2500ms, and outside of MAD criteria
    if ~isfolder ([oPath filesep '1_post_dbl_RT'])
        mkdir ([oPath filesep '1_post_dbl_RT']);
    end

    % Create folder for concatenated EEG files
    if ~isfolder ([oPath filesep '2_preprocessed'])
        mkdir ([oPath filesep '2_preprocessed']);
    end

    % Create folder for analysis-ready response-locked data with all participants
    if ~isfolder ([oPath filesep '3_processed_response'])
        mkdir ([oPath filesep '3_processed_response']);
    end

    % Create folder for analysis-ready stimulus-locked data with all participants
    % This is not used for the first manuscript, but is done to prepare for
    % stimulus-locked EEG analyses
    if ~isfolder ([oPath filesep '4_processed_stimulus'])
        mkdir ([oPath filesep '4_processed_stimulus']);
    end

    % Create folder to track preprocessing behavioral outputs
    if ~isfolder ([oPath filesep 'preprocess_1_tracking'])
        mkdir ([oPath filesep 'preprocess_1_tracking']);
    end

    % Create folder to keep behavioral info after removal of double-clicking
    if ~isfolder ([oPath filesep 'preprocess_2_behavioral_dbl'])
        mkdir ([oPath filesep 'preprocess_2_behavioral_dbl']);
    end

    % Create folder to keep behavioral info after removal of double-clicking
    % and basd on RT threshold
    if ~isfolder ([oPath filesep 'preprocess_3_behavioral_postRT'])
        mkdir ([oPath filesep 'preprocess_3_behavioral_postRT']);
    end

    % Initialize overall metric tracking
    File_length = [];
    Number_chan_kept = [];
    Percent_chan_kept=[];
    Name_chan_removed = [];
    Number_ICA_removed = [];
    Percent_ICA_removed = [];
    Num_trials_rejected_perRT = [];
    Num_trials_kept_intermed = [];
    Num_trials_kept_final= [];
    Num_errors_kept_final = [];
    Accuracy_dbl = [];
    Correct_dbl = [];
    Error_dbl = [];


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%% Begin preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Processing each participant separately
    for i = 1:length(Participants)
        i/length(Participants)

        % Initialize EEG objects, just in case so there is no "re-using" of
        % previous participant
        EEG_beh = [];
        EEG_pre = [];
        EEG_resp = [];
        EEG_resp2 = [];
        EEG_stim = [];
        EEG_stim2 = [];

        eeglab

        % Set random number generator seed for reproduceability
        rng(0, 'twister')

        % Load file and make sampling rate 500 Hz
        EEG_beh = pop_loadbv(iPath, Participants{i});
        Name = Participants{i};
        Name = Name(1:end-5);
        if EEG_beh.srate > 500
            EEG_beh = pop_resample( EEG_beh, 500);
        end

        EEG_beh = eeg_checkset(EEG_beh);
        File_length(i) = EEG_beh.xmax;

        %%%%%%%%%%% Process behavioral data first %%%%%%%%%%%%%%%%%
        %% Remove breaks, defined as 5 seconds of empty no event before and after an event
        % 5 seconds, becuase the epochs for analyses are about -3 to 3.5
        % seconds, so giving sufficient buffer before and after to avoid
        % boundary events interfering

        % Adapted from https://sccn.ucsd.edu/pipermail/eeglablist/2014/007596.html
        % By James Jones-Rounds
        indices_to_delete = [];
        break_in_sec = 5;
        break_in_pts = break_in_sec * EEG_beh.srate;

        % Loop through each event
        for ev_i =  1:length(EEG_beh.event)
            prev_ev = [];
            current_ev = [];
            next_ev = [];

            % Find event latencies in points
            % For first event
            if ev_i == 1
                % prev_ev = [];
                current_ev = EEG_beh.event(ev_i).latency;
                next_ev = EEG_beh.event(ev_i + 1).latency;

                % calculate latency from previous and next events
                prev_inter_trial_latency = current_ev - 1;
                next_inter_trial_latency = next_ev - current_ev;

                % track points to delete
                % This takes care of any 5+ second long break before the
                % task begins
                if  prev_inter_trial_latency > break_in_pts
                    indices_to_delete = [indices_to_delete; 1 (current_ev - break_in_pts)];
                end
                % This takes care of any 10+ second long break between 1st
                % and 2nd event
                if next_inter_trial_latency > break_in_pts*2
                    indices_to_delete = [indices_to_delete; (current_ev + break_in_pts) (next_ev - break_in_pts)];
                end


                % For last event
            elseif ev_i == length(EEG_beh.event)
                % This takes care of any 5+ second break after the final
                % trial

                % calculate latency from previous event
                % prev_ev = EEG_beh.event(ev_i - 1).latency;
                current_ev = EEG_beh.event(ev_i).latency;
                next_ev = [];

                % prev_inter_trial_latency = current_ev - prev_ev;
                next_inter_trial_latency = length(EEG_beh.times) - current_ev;

                % if  prev_inter_trial_latency > break_in_pts*2
                %     indices_to_delete = [indices_to_delete; (prev_ev + break_in_pts) (current_ev - break_in_pts)];
                % end
                if next_inter_trial_latency > break_in_pts
                    indices_to_delete = [indices_to_delete; (current_ev + break_in_pts) length(EEG_beh.times)];
                end


                % For all other events
                % This takes care of any 10+ second breaks between this and
                % the next trial
            else
                % prev_ev = EEG_beh.event(ev_i - 1).latency;
                current_ev = EEG_beh.event(ev_i).latency;
                next_ev = EEG_beh.event(ev_i + 1).latency;

                % calculate latency from previous and next events
                % prev_inter_trial_latency = current_ev - prev_ev;
                next_inter_trial_latency = next_ev - current_ev;

                % track points to delete
                % if  prev_inter_trial_latency > break_in_pts*2
                %     indices_to_delete = [indices_to_delete; (prev_ev + break_in_pts) (current_ev - break_in_pts)];
                % end
                if next_inter_trial_latency > break_in_pts*2
                    indices_to_delete = [indices_to_delete; (current_ev + break_in_pts) (next_ev - break_in_pts)];
                end
            end
        end

        % actually delete events, and only if there is something to
        % delete
        if length(indices_to_delete) > 0
            EEG_beh = pop_select(EEG_beh, 'nopoint', indices_to_delete);
        end

        %% Remove trials based on responses

        % Remove second of double-responses
        [Resp_dbl,PStimUR_dbl,PStim_dbl,RespLat_dbl] = eeg_context(EEG_beh,Resp_con,Stim_con,-1);  % find latency from previous "stimulus" event
        [Double_dbl,PRespUR_dbl,PResp_dbl,DoublLat_dbl] = eeg_context(EEG_beh,Resp_con,Resp_con,-1);  % find latency from previous "resposne" event
        rejev_dbl = []; % rejected events

        for e=1:length(Resp_dbl) % For each response,
            if (DoublLat_dbl(e) < RespLat_dbl(e) | isnan(DoublLat_dbl(e)))
                % Remove intervening 'Response' event (included only if latency from
                % previous "response" is longer than latency from previous "stimulus")
                % Equality in opposite direction, becuase info exported as
                % negative values
                % isnan necessary for first event
            else
                rejev_dbl = [rejev_dbl e];  % otherwise reject (these are 2nd of double responses
            end
        end

        rejev_dbl_ev = Resp_dbl(rejev_dbl); %select the event numbers to reject
        EEG_beh = pop_selectevent( EEG_beh, 'omitevent',rejev_dbl_ev ,'deleteevents','on'); %actually reject
        events_e = eeg_eventtable(EEG_beh, 'unit', 'seconds', 'exportFile', [oPath filesep 'preprocess_2_behavioral_dbl' filesep Name '_dbl.csv']);  % save behavioral data

        % Calculate accuracy
        Error_dbl(i) = length(find(ismember({EEG_beh.event.type}, Resp_Err))); % find position of error responses
        Correct_dbl(i) = length(find(ismember({EEG_beh.event.type}, Resp_Cor))); % find position of correct responses
        Accuracy_dbl(i) = (Correct_dbl(i) / (Error_dbl(i) + Correct_dbl(i)));
        EEG_beh.etc.correct_dbl = Correct_dbl(i);
        EEG_beh.etc.error_dbl = Error_dbl(i);
        EEG_beh.etc.accuracy = Accuracy_dbl(i);

        %This removes events based on RT
        [Resp_rt,PStimUR_rt,PStim_rt,RespLat_rt] = eeg_context(EEG_beh,Resp_con,Stim_con,-1);  % find latency from previous "stimulus" event
        rejev_rt = [];

        % Determine threshold for noise (median absolute deviation; equivalent
        % to 3 standard deviations). This is to remove any trials that
        % may not be assessing what we are invetigating (e.g.,
        % reflex, inattention, spacing). Ideally, MAD would take care
        % of any within-participant outliers.
        % However, added absolute values to ensure
        % Too quick response defined as: below MAD criteria OR 100ms
        % Too slow response defined as: above MAD criteria OR 1500ms
        znoise = mad(RespLat_rt,1)*1.4826;
        high_MAD = median(RespLat_rt, 'omitnan') - znoise * 3;
        low_MAD = median(RespLat_rt,'omitnan')  + znoise * 3;

        % Based on observation that most papers
        high_threshold = max(high_MAD, -1500); % 1500ms or faster RT
        low_threshold = min(low_MAD, -100); % 100ms or slower RT

        for f=1:length(RespLat_rt) % For each response,
            if RespLat_rt(f) > high_threshold & RespLat_rt(f) < low_threshold
                % if RT in acceptable range (negative because calculated
                % as [stimuli - response])
            else
                rejev_rt = [rejev_rt f];  % otherwise reject (these are 2nd of double responses)
            end
        end

        Num_trials_rejected_perRT(i) = length(rejev_rt);
        rejev_rt_ev = Resp_rt(rejev_rt); % select the event numbers to reject
        EEG_beh = pop_selectevent( EEG_beh, 'omitevent',rejev_rt_ev ,'deleteevents','on'); % actually reject
        events_f = eeg_eventtable(EEG_beh, 'unit', 'seconds', 'exportFile', [oPath filesep 'preprocess_3_behavioral_postRT' filesep Name '_rt.csv']); % save behavioral data

        % Save EEG file cleaned for double-responses and RT
        pop_saveset(EEG_beh, [Name '_post_dbl_RT'], [oPath filesep '1_post_dbl_RT']);


        %% %%%%%%%%%%% Preprocess EEG %%%%%%%%%%%
        EEG_pre = pop_select(EEG_beh,'nochannel',{'VEOGL'});

        % Add online reference Cz
        EEG_pre.data(end+1,:) = 0;
        EEG_pre.nbchan = size(EEG_pre.data,1);
        if ~isempty(EEG_pre.chanlocs)
            EEG_pre.chanlocs(end+1).labels = 'Cz';
        end
        [ALLEEG EEG_pre CURRENTSET] = eeg_store(ALLEEG, EEG_pre, CURRENTSET);
        EEG_pre = pop_chanedit( EEG_pre,'lookup',chan_locations,'setref',{'1:63','Cz'});

        % Offine re-reference to average
        EEG_pre = pop_reref( EEG_pre, []);
        EEG_pre = eeg_checkset( EEG_pre );

        EEG_pre.full_chan = EEG_pre.chanlocs;
        original_channel_labels = {EEG_pre.full_chan.labels};
        EEG_pre.chaninfo.removedchans = [];

        % Filter the data with 0.1hz highpass and 249hz lowpass
        EEG_pre = pop_eegfiltnew(EEG_pre, 0.1, 249,[],0,[],0);
        EEG_pre = eeg_checkset( EEG_pre );

        % Reduce line noise in the data
        EEG_pre = pop_cleanline(EEG_pre, 'bandwidth',2,'chanlist',[1:EEG_pre.nbchan],'computepower',1,'linefreqs',...
            [60 120 180 240] ,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype',...
            'Channels','tau',100,'verb',0,'winsize',4,'winstep',1, 'ComputeSpectralPower','False');
        EEG_pre = eeg_checkset(EEG_pre);

        % Remove bad channels defined by 5 sec flatline and
        % below minimum correlation with nearby electrodes
        EEG_pre = pop_clean_rawdata(EEG_pre, 'FlatlineCriterion', 4,'ChannelCriterion', Threshold_Min_Correlation,'LineNoiseCriterion','off','Highpass','off', ...
            'BurstCriterion','off','WindowCriterion','off','BurstRejection','off','Distance','Euclidian');
        selected_channel_locations=EEG_pre.chanlocs;

        % Save the names of the rejected channels for output table after the pipeline finishes
        selected_channel_labels={selected_channel_locations.labels};
        bad_channels_removed= setdiff(original_channel_labels, selected_channel_labels);


        % Rereference to average (with interpolation of channels, and removal after referencing)
        if length(bad_channels_removed) > 0
            EEG_pre = pop_reref( EEG_pre,[],'interpchan',[]);
        else
            EEG_pre = pop_reref( EEG_pre,[]);
        end

        % Save outputs
        Number_chan_original(i) = size(original_channel_labels,2);
        Number_chan_kept(i) = size(selected_channel_locations,2);

        if ~isempty(bad_channels_removed)
            Name_chan_removed{i} = [sprintf('%s ',bad_channels_removed{1:end-1}),bad_channels_removed{end}];
        else
            Name_chan_removed{i} = [];
        end
        Percent_chan_kept(i)=Number_chan_kept(i)/Number_chan_original(i)* 100;

        % Remove bad segments, using EEGLAB ARS
        EEG_pre = pop_clean_rawdata(EEG_pre, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off', ...
            'BurstCriterion',Threshold_ARS,'WindowCriterion','off','BurstRejection','on','Distance','Euclidian');
        % Separated Bad channel detection and bad segment removals to add a
        % step to re-reference with interpolation to reduce effect of bad
        % channel.
        % 'WindowCriterionTolerances',[-Inf 7] removed. These were too
        % sensitive and removed many blinks even with 100.
        % 'BurstCriterion', 20 changed to 100. 20 was too sensitive and
        % removed many blinks.

        % Run ICA, IC Label, and reject components
        EEG_pre = pop_runica(EEG_pre, 'extended',1);
        EEG_pre = pop_iclabel(EEG_pre, 'default');
        EEG_pre = pop_icflag(EEG_pre, [NaN NaN;Threshold_ICA 1;Threshold_ICA 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]);
        % Thresholds changed from .90 to .75, due to some blinks being picked
        % up as multiple components, and each labeled between .75 and .90. .80
        % was tried, but there was at least 1 participant with component 1 .793
        % blink (with next being a pretty big muscle probability) and
        % non-removal led to deletion of 90%+ trials being rejeted.
        Number_ICA_removed(i) = sum(EEG_pre.reject.gcompreject);
        Percent_ICA_removed(i)=Number_ICA_removed(i)/Number_chan_kept(i)* 100;
        EEG_pre = pop_subcomp( EEG_pre, [], 0);

        % Interpolate removed channels
        EEG_pre = pop_interp(EEG_pre, EEG_pre.full_chan);

        % Save Data
        EEG_pre = pop_saveset(EEG_pre, 'filename',[Name,'_preprocessed.set'],'filepath',[oPath filesep '2_preprocessed']);

        %% %%%%%%% Preparation for TF/ERP analyses %%%%%%
        % Create files for response-locked analyses
        EEG_resp = pop_loadset([Name,'_preprocessed.set'], [oPath filesep '2_preprocessed']);
        EEG_resp = pop_epoch(EEG_resp, Resp_con, epoch_window);
        Num_trials_kept_intermed(i) = length(EEG_resp.epoch);

        EEG_resp = pop_rmbase(EEG_resp, [-200 -50]);   % BL for response locked ERPs/absolute voltage
        EEG_resp = pop_selectevent( EEG_resp, 'latency','-3000 <= 1','type', {'S  1', 'S  2', 'S  3', 'S  4', 'S  9', 'S 10'},...
            'deleteevents','on','deleteepochs','off','invertepochs','off');
        [EEG_resp Indexes] = pop_eegthresh( EEG_resp, 1, [1:EEG_resp.nbchan], -Threshold_Trials, Threshold_Trials, -0, 0.5, 0, 1);
        good_events = [EEG_resp.event(:).urevent];
        Num_trials_kept_final(i) = length(EEG_resp.epoch);
        Num_errors_kept_final(i) = length(find(ismember({EEG_resp.event.type}, Resp_Err)));

        EEG_resp2 = pop_loadset([Name,'_preprocessed.set'], [oPath filesep '2_preprocessed']);
        idx_resp = setdiff([EEG_resp2.event.urevent], good_events);
        for h = 1:length(EEG_resp2.event)
            if  ismember(EEG_resp2.event(h).type, Resp_con) && ismember(EEG_resp2.event(h).urevent, idx_resp)
                EEG_resp2.event(h).type = '98';
            end
        end

        pop_saveset(EEG_resp2, [Name '_processed_response'], [oPath filesep '3_processed_response']); % Save the file as a new EEGLAB dataset


        % Create files for stimulus-locked analyses
        eeglab
        EEG_stim = pop_loadset([Name,'_preprocessed.set'], [oPath filesep '2_preprocessed']);
        EEG_stim = pop_epoch( EEG_stim, Stim_con, epoch_window);

        EEG_stim = pop_rmbase( EEG_stim, [-200 -50]);   % BL for stimulus locked ERPs/absolute voltage
        EEG_stim = pop_selectevent( EEG_stim, 'latency','-1 <= 3000','type', {'S  1', 'S  2', 'S  3', 'S  4', 'S  9', 'S 10'},...
            'deleteevents','on','deleteepochs','off','invertepochs','off');
        [EEG_stim Indexes] = pop_eegthresh( EEG_stim, 1, [1:EEG_stim.nbchan], -Threshold_Trials, Threshold_Trials, -0, 0.5, 0, 1);
        good_events = [EEG_stim.event(:).urevent];

        EEG_stim2 = pop_loadset([Name,'_preprocessed.set'], [oPath filesep '2_preprocessed']);
        idx_stim = setdiff([EEG_stim2.event.urevent], good_events);
        for g = 1:length(EEG_stim2.event)
            if  ismember(EEG_stim2.event(g).type, Stim_con) && ismember(EEG_stim2.event(g).urevent, idx_stim)
                EEG_stim2.event(g).type = '98';
            end
        end

        pop_saveset(EEG_stim2, [Name '_processed_stimulus'], [oPath filesep '4_processed_stimulus']); % Save the file as a new EEGLAB dataset

        table_all = table({Participants{1:i}}',File_length',Number_chan_kept',Percent_chan_kept',...
            Name_chan_removed', Number_ICA_removed', Percent_ICA_removed',...
            Num_trials_rejected_perRT', Num_trials_kept_intermed', Num_trials_kept_final', Num_errors_kept_final',...
            Accuracy_dbl', Correct_dbl', Error_dbl');

        table_all.Properties.VariableNames ={'Participants', 'File_length','Number_chan_kept','Percent_chan_kept',...
            'Name_chan_removed', 'Number_ICA_removed', 'Percent_ICA_removed',...
            'Num_trials_rejected_perRT','Num_trials_kept_intermed','Num_trials_kept_final', 'Num_errors_kept_final',...
            'Accuracy', 'Correct_num', 'Error_num'};

        cd([oPath filesep 'preprocess_1_tracking'])
        writetable(table_all, [oPath filesep 'preprocess_1_tracking' filesep 'Preprocessing_trials_' datestr(now,'yyyy-mm-dd'),'.csv']);

        toc
    end
end
fprintf('Step 1 Preprocessing is done \n');
