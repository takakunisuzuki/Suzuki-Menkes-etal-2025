%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Response Monitoring Theta-Band Activities across Emotional Contexts in 
% Schizophrenia- and Bipolar-Spectrum Disorders
% Suzuki, Menkes, et al.
%
% Script to conduct Steps 2-4 using an earlier baseline (-200 to -400ms) for ERPs. 
%
% Completed using MATLAB 2024a & EEGLAB 2024.0, on Windows 10 Enterprise
%
% Author: Takakuni Suzuki
% First drafted June 2023
% Last updated January 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

%%%%%%%%% Setting up the preprocessing parameters, folders, etc. %%%%%%%%%
ParentPath = ''; % Specifcy data folder
folders = {'FARR' 'FNEG' 'FPOS'};

%% Step 2 equivalent
for fol = 1:length(folders)
    % Set paths
    TaskPath = [ParentPath filesep folders{fol}];
    iPath = [TaskPath filesep 'Preprocess_Auto' filesep '3_processed_response'];
    oParentPath = [TaskPath filesep 'Output_Auto_updated'];
    o1Path = [oParentPath filesep 'Response_Matched_EBL'];
    o2Path = [oParentPath filesep 'Response_PowerITPC_EBL'];

    % Set suffix
    isuffix = '_processed_response.set';
    osuffix = '_matched_response_EBL.set';

    % Enter response marker info
    Resp_Cor = {'S  3' 'S  4'};
    Resp_Err = {'S  9' 'S 10'};
    Resp_Events = [Resp_Cor, Resp_Err];
    Stim_Events = {'S  1' 'S  2'};

    window = [-3 (-3 + (500 * 5.5 - 1)/500)];

    %% Create folders
    if ~isdir ([oParentPath])
        mkdir ([oParentPath]);
    end

    if ~isdir ([o1Path])
        mkdir ([o1Path]);
    end

    if ~isdir ([o2Path])
        mkdir ([o2Path]);
    end

    %% %%%%% Start extraction %%%%%
    % Read in participant files
    cd(iPath);
    Participants=dir(['*' isuffix]);
    Participants={Participants.name};

    % Remove extra parts of the file names
    Participants_n = string(nan(length(Participants), 1));
    for j = 1:length(Participants)
        name = char(Participants(j));
        Participants_n(j) = string(name(1:13));
    end

    % Start with matching number of correct trails to number of error trials
    for i = 1:length(Participants_n)
        % Initialize data matrix
        data = [];
        % row = [];
        EEG = [];
        % EEG_base = [];
        EEG_correct= [];
        EEG_correct_non = [];
        EEG_error = [];

        cd(iPath);
        eeglab
        EEG = pop_loadset([Participants_n{i}, isuffix], iPath);

        % Create a new structure to claculate reaction times
        % (only taking differences here)
        for k = 2:length(EEG.event)
            EEG.event(k).rt = EEG.event(k).latency - EEG.event(k-1).latency;
        end

        % Epoch
        EEG = pop_epoch( EEG, Resp_Events, window);

        % To remvoe event codes from previous and after the responses
        EEG = pop_selectevent( EEG, 'latency','-100 <= 1','type', {'S  3', 'S  4', 'S  9', 'S 10'},...
            'deleteevents','on','deleteepochs','off','invertepochs','off');

        % Extract indices for error and correct trials
        idx_error = find(ismember({EEG.event.type}, Resp_Err));
        idx_correct = find(ismember({EEG.event.type}, Resp_Cor));

        %Initialize time frequency data
        % basepower = NaN(EEG.nbchan,num_frex);
        % eegpower_error = NaN(EEG.nbchan,num_frex,EEG.pnts);
        % eegpower_correct= NaN(EEG.nbchan,num_frex,EEG.pnts);
        % eegitpc_error = NaN(EEG.nbchan,num_frex,EEG.pnts);
        % eegitpc_correct = NaN(EEG.nbchan,num_frex,EEG.pnts);
        erp_error = NaN((EEG.nbchan-2),EEG.pnts);
        erp_correct = NaN((EEG.nbchan-2),EEG.pnts);
        erp_correct_non = NaN((EEG.nbchan-2),EEG.pnts);


        if length(idx_correct) > 0 && length(idx_error) > 0 && EEG.etc.accuracy > .50

            % Initialize final trial matrix
            idxAll = zeros(length(idx_error), 2);
            idxAll(:,1) = idx_error;

            % Start loop to extract the indices to keep
            % Loop through all error trials
            for e = 1:length(idx_error)
                % Something large to compare the aboslute differnces in RT
                diff = 9999;
                % Start going through each correct trials...
                for c = 1:length(idx_correct)
                    % Except those already used
                    if ~ismember(idx_correct(c), idxAll(:,2))
                        % Find absolute difference and compare to previous largest
                        % absolute difference
                        if abs(EEG.event(idx_error(e)).rt - EEG.event(idx_correct(c)).rt) < diff
                            diff = abs(EEG.event(idx_error(e)).rt - EEG.event(idx_correct(c)).rt);
                            idxAll(e,2) = idx_correct(c);
                        end
                    end
                end
            end

            % Initialize epoch keeping
            idxEpo = zeros(1, length(idxAll(:)));
            idxEpo_base = zeros(1, length(idxAll(:)));
            for p = 1:length(idxAll(:))
                idxEpo(p) = EEG.event(idxAll(p)).epoch;
                % Saving corresponding stimulus for each epoch for later
                % baseline extraction
                idxEpo_base(p) = (EEG.event(idxAll(p)).urevent-1);
            end

            % Reject unneeded epochs
            EEG = pop_selectevent( EEG, 'epoch', idxEpo ,'deleteevents','off','deleteepochs','on','invertepochs','off');

            % Save the data that only keeps correct trials matched to error trials.
            pop_saveset(EEG, [Participants_n{i} osuffix], o1Path);

            %% Start time-frequency decomposition
            %% TF of baseline
            % EEG_base = pop_loadset([Participants_n{i}, isuffix], iPath);
            %
            % idx_base = setdiff([EEG_base.event.urevent], idxEpo_base);
            % for h = 1:length(EEG_base.event)
            %     if  ismember(EEG_base.event(h).type, Stim_Events) && ismember(EEG_base.event(h).urevent, idx_base)
            %         EEG_base.event(h).type = '98';
            %     end
            % end
            %
            % % Extract epochs based on stimulus presentation
            % % Extracting 300ms extra on the front end for baseline window
            % EEG_base = pop_epoch( EEG_base, Stim_Events, [-1.8 1.51]);
            % EEG_base = eeg_checkset( EEG_base );
            %
            % % Define wavelet parameters
            % time = -1:1/EEG_base.srate:1;
            % frex = logspace(log10(min_freq),log10(max_freq),num_frex);
            % s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
            %
            % % Define convolution parameters
            % n_wavelet_base            = length(time);
            % n_data_base               = EEG_base.pnts*EEG_base.trials;
            % n_convolution_base        = n_wavelet_base+n_data_base-1;
            % n_conv_pow2_base          = pow2(nextpow2(n_convolution_base));
            % half_of_wavelet_size_base = (n_wavelet_base-1)/2;
            %
            % % Loop through each electrode
            % for k = 1:EEG_base.nbchan
            %     chan2use_base = EEG_base.chanlocs(k).labels;
            %
            %     % get FFT of data
            %     eegfft_base = fft(reshape(EEG_base.data(strcmpi(chan2use_base,{EEG_base.chanlocs.labels}),:,:),1,EEG_base.pnts*EEG_base.trials),n_conv_pow2_base);
            %
            %     % Find locations of the baseline window rangea
            %     baseidx = dsearchn(EEG_base.times',[-300 -50]');
            %
            %     % loop through frequencies and compute synchronization
            %     for fi=1:num_frex
            %         wavelet_base = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2_base );
            %
            %         % convolution
            %         eegconv_base = ifft(wavelet_base.*eegfft_base);
            %         eegconv_base = eegconv_base(1:n_convolution_base);
            %         eegconv_base = eegconv_base(half_of_wavelet_size_base+1:end-half_of_wavelet_size_base);
            %
            %         % Average power over trials
            %         temppower_base = mean(abs(reshape(eegconv_base,EEG_base.pnts,EEG_base.trials)).^2,2);
            %         basepower(k,fi,:) = mean(temppower_base(baseidx(1):baseidx(2)));
            %     end
            % end

            %% TF for error trials
            EEG_error = pop_loadset([Participants_n{i} osuffix], o1Path);

            % extract only error trials
            EEG_error = pop_selectevent( EEG_error, 'type',Resp_Err,'deleteevents','off','deleteepochs','on','invertepochs','off');
            EEG_error = eeg_checkset( EEG_error );

            % % define wavelet parameters
            % time = -1:1/EEG_error.srate:1;
            % frex = logspace(log10(min_freq),log10(max_freq),num_frex);
            % s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
            %
            % % definte convolution parameters
            % n_wavelet_error            = length(time);
            % n_data_error               = EEG_error.pnts*EEG_error.trials;
            % n_convolution_error        = n_wavelet_error+n_data_error-1;
            % n_conv_pow2_error          = pow2(nextpow2(n_convolution_error));
            % half_of_wavelet_size_error = (n_wavelet_error-1)/2;
            %
            % % Loop through each electrode
            % for j = 1:EEG_error.nbchan %EEG.nbchan
            %     chan2use_error = EEG_error.chanlocs(j).labels;
            %
            %     % get FFT of data
            %     eegfft_error = fft(reshape(EEG_error.data(strcmpi(chan2use_error,{EEG_error.chanlocs.labels}),:,:),1,EEG_error.pnts*EEG_error.trials),n_conv_pow2_error);
            %
            %     % loop through frequencies and compute synchronization
            %     for fi=1:num_frex
            %         wavelet_error = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2_error );
            %
            %         % convolution
            %         eegconv_error = ifft(wavelet_error.*eegfft_error);
            %         eegconv_error = eegconv_error(1:n_convolution_error);
            %         eegconv_error = eegconv_error(half_of_wavelet_size_error+1:end-half_of_wavelet_size_error);
            %
            %         % Average power over trials and adjust for baseline
            %         temppower_error = mean(abs(reshape(eegconv_error, EEG_error.pnts, EEG_error.trials)).^2,2);
            %         eegpower_error(j,fi,:) = 10*log10(temppower_error./basepower(j,fi));
            %         eegitpc_error(j,fi,:) = abs(mean(exp(1i*angle(reshape(eegconv_error,EEG_error.pnts,EEG_error.trials))),2));
            %     end
            % end
            % tf_chan = EEG_error.chanlocs;


            % For ERN
            EEG_error = pop_reref( EEG_error, [10 21] );% Mastoid (tp9/tp10) reference
            EEG_error = pop_rmbase( EEG_error, [-400 -200]); %Baseline
            for k = 1:EEG_error.nbchan %EEG.nbchan
                erp_error(k,:) = nanmean(EEG_error.data(k,:,:),3) ;
            end
            erp_chan = EEG_error.chanlocs;

            %% TF for Correct trials
            EEG_correct = pop_loadset([Participants_n{i} osuffix], o1Path);

            % % define wavelet parameters
            % time = -1:1/EEG_correct.srate:1;
            % frex = logspace(log10(min_freq),log10(max_freq),num_frex);
            % s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
            %
            % % extract only correct trials
            % EEG_correct = pop_selectevent( EEG_correct, 'type',Resp_Cor,'deleteevents','off','deleteepochs','on','invertepochs','off');
            % EEG_correct = eeg_checkset( EEG_correct );
            %
            % % Define convolution parameters
            % n_wavelet_correct            = length(time);
            % n_data_correct               = EEG_correct.pnts*EEG_correct.trials;
            % n_convolution_correct        = n_wavelet_correct+n_data_correct-1;
            % n_conv_pow2_correct          = pow2(nextpow2(n_convolution_correct));
            % half_of_wavelet_size_correct = (n_wavelet_correct-1)/2;
            %
            % % Loop through each electrode
            % for l = 1:EEG_correct.nbchan %EEG.nbchan
            %     chan2use_correct = EEG_correct.chanlocs(l).labels;
            %
            %     % get FFT of data
            %     eegfft_correct = fft(reshape(EEG_correct.data(strcmpi(chan2use_correct,{EEG_correct.chanlocs.labels}),:,:),1,EEG_correct.pnts*EEG_correct.trials),n_conv_pow2_correct);
            %
            %     % loop through frequencies and compute synchronization
            %     for fi=1:num_frex
            %         wavelet_correct = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2_correct );
            %
            %         % convolution
            %         eegconv_correct = ifft(wavelet_correct.*eegfft_correct);
            %         eegconv_correct = eegconv_correct(1:n_convolution_correct);
            %         eegconv_correct = eegconv_correct(half_of_wavelet_size_correct+1:end-half_of_wavelet_size_correct);
            %
            %         % Average power over trials and adjust for baseline
            %         temppower_correct = mean(abs(reshape(eegconv_correct,EEG_correct.pnts,EEG_correct.trials)).^2,2);
            %         eegpower_correct(l,fi,:) = 10*log10(temppower_correct./basepower(l,fi));
            %         eegitpc_correct(l,fi,:) = abs(mean(exp(1i*angle(reshape(eegconv_correct,EEG_correct.pnts,EEG_correct.trials))),2));
            %     end
            % end

            % For matched CRN
            EEG_correct = pop_reref( EEG_correct, [10 21] );% Mastoid (tp9/tp10) reference
            EEG_correct = pop_rmbase( EEG_correct, [-400 -200]); %Baseline
            for m = 1:EEG_correct.nbchan %EEG.nbchan
                erp_correct(m,:) = nanmean(EEG_correct.data(m,:,:),3) ;
            end

            % For non-matched CRN
            EEG_correct_non = pop_loadset([Participants_n{i}, isuffix], iPath);
            EEG_correct_non = pop_epoch( EEG_correct_non, Resp_Cor, window);
            EEG_correct_non = pop_reref( EEG_correct_non, [10 21] );% Mastoid (tp9/tp10) reference
            EEG_correct_non = pop_rmbase( EEG_correct_non, [-400 -200]); %Baseline
            for n = 1:EEG_correct_non.nbchan
                erp_correct_non(n,:) = nanmean(EEG_correct_non.data(n,:,:),3) ;
            end
        end

        % data.power.error = eegpower_error;
        % data.power.correct = eegpower_correct;
        % data.ITPC.error = eegitpc_error;
        % data.ITPC.correct = eegitpc_correct;
        data.erp.error = erp_error;
        data.erp.correct = erp_correct;
        data.erp.correct_non = erp_correct_non;
        % data.freq = frex;
        data.time = EEG.times;
        % data.chan_tf = tf_chan;
        data.chan_erp = erp_chan;
        data.UsableError = length(idx_error);
        data.UsableCorrect = length(idx_correct);
        data.AllError = EEG.etc.error_dbl;
        data.AllCorrect = EEG.etc.correct_dbl;
        data.accuracy = EEG.etc.accuracy;

        %Attaching group info
        if strcmp(Participants_n{i}(1:10),strcat('TNA_', folders{fol}, '_1'))
            data.group = 'HC';
        else
            if strcmp(Participants_n{i}(1:10),strcat('TNA_', folders{fol}, '_2'))
                data.group = 'SZ';
            else
                if strcmp(Participants_n{i}(1:10),strcat('TNA_', folders{fol}, '_3'))
                    data.group = 'BD';
                else
                    data.group = [];
                end
            end
        end

        save(char(strcat(o2Path, filesep, Participants_n(i), '_data_EBL', '.mat')), 'data');
    end
end


%% Step 3 equivalent
for fol = 1:length(folders)
    TaskPath = [ParentPath filesep folders{fol} filesep 'Output_Auto_updated'];

    %% Define parameters
    gr = {'HC','SZ','BD'};

    %% Define paths and create folders
    iPath = [TaskPath filesep 'Response_PowerITPC_EBL'];
    oPath = [TaskPath filesep 'AnalysisData_EBL'];

    if ~isdir ([oPath])
        mkdir ([oPath]);
    end


    %% Concatenate TF and ERP data
    % Load files with "_data.mat" suffix in the
    cd(iPath);
    files = dir('*_data_EBL.mat');
    filenames = {files.name};
    Participants = regexprep(filenames, '_data_EBL.mat', '');

    % Import one dataset to extract parameters (e.g., channel locations)
    load(char(strcat([Participants(1)], '_data_EBL.mat')))
    % num_chan_tf = length(data.chan_tf);
    num_chan_erp = length(data.chan_erp);
    % num_frex = length(data.freq);
    num_time = length(data.time);

    % Initiate TF and ERP data matrices for all participants
    % Power_E_All(:,:,:) = nan(num_chan_tf, num_frex, num_time);
    % Power_C_All(:,:,:) = nan(num_chan_tf, num_frex, num_time);
    % ITPC_E_All(:,:,:) =  nan(num_chan_tf, num_frex, num_time);
    % ITPC_C_All(:,:,:) =  nan(num_chan_tf, num_frex, num_time);
    ERP_E_All(:,:) = nan(num_chan_erp, num_time);
    ERP_C_All(:,:) = nan(num_chan_erp, num_time);
    ERP_C_non_All(:,:) = nan(num_chan_erp, num_time);

    % Import TF analysis and ERP data from all participants
    % going through channels to avoid memory issues
    for l = 1:num_chan_erp % loop through all channels
        % Initialize matrix for this channel's data
        % Power_E_All_temp = nan(length(Participants), num_frex, num_time);
        % Power_C_All_temp = nan(length(Participants), num_frex, num_time);
        % ITPC_E_All_temp = nan(length(Participants), num_frex, num_time);
        % ITPC_C_All_temp = nan(length(Participants), num_frex, num_time);
        ERP_E_All_temp = nan(length(Participants), num_time);
        ERP_C_All_temp = nan(length(Participants), num_time);
        ERP_C_All_non_temp = nan(length(Participants),  num_time);

        for i = 1:length(Participants) % loop through all participants for each channel
            load(char(strcat([Participants(i)], '_data_EBL.mat')))

            % Use data with more 8 or more usable errors
            if data.UsableError > 7
                % Power_E_All_temp(i,:,:) = data.power.error(l,:,:);
                % Power_C_All_temp(i,:,:) = data.power.correct(l,:,:);
                % ITPC_E_All_temp(i,:,:) = data.ITPC.error(l,:,:);
                % ITPC_C_All_temp(i,:,:) = data.ITPC.correct(l,:,:);

                % if l <= num_chan_erp % ERP channels have 2 less channels due to referencing to mastoid-equivalents
                ERP_E_All_temp(i,:) = data.erp.error(l,:);
                ERP_C_All_temp(i,:) = data.erp.correct(l,:);
                ERP_C_All_non_temp(i,:) = data.erp.correct_non(l,:);
                % End % end ERP 
            end % End usable data extraction
        end % End participant loop for this particular channel

        % Calculate mean across participants for this channel
        % Power_E_All(l,:,:) = squeeze(mean(Power_E_All_temp(:,:,:), 1, 'omitnan'));
        % Power_C_All(l,:,:) = squeeze(mean(Power_C_All_temp(:,:,:), 1, 'omitnan'));
        % ITPC_E_All(l,:,:) =  squeeze(mean(ITPC_E_All_temp(:,:,:), 1, 'omitnan'));
        % ITPC_C_All(l,:,:) =  squeeze(mean(ITPC_C_All_temp(:,:,:), 1, 'omitnan'));

        % if l <= num_chan_erp
        ERP_E_All(l,:) = squeeze(mean(ERP_E_All_temp(:,:), 1, 'omitnan'));
        ERP_C_All(l,:)  = squeeze(mean(ERP_C_All_temp(:,:), 1, 'omitnan'));
        ERP_C_non_All(l,:)  = squeeze(mean(ERP_C_All_non_temp(:,:), 1, 'omitnan'));
        % end
    end

    % Save data with average from all participants
    % save(char(strcat(oPath, filesep, 'Power_E_all.mat')), 'Power_E_All')
    % save(char(strcat(oPath, filesep, 'Power_C_all.mat')), 'Power_C_All')
    % save(char(strcat(oPath, filesep, 'ITPC_E_all.mat')), 'ITPC_E_All')
    % save(char(strcat(oPath, filesep, 'ITPC_C_all.mat')), 'ITPC_C_All')
    save(char(strcat(oPath, filesep, 'ERP_E_EBL_all.mat')), 'ERP_E_All')
    save(char(strcat(oPath, filesep, 'ERP_C_EBL_all.mat')), 'ERP_C_All')
    save(char(strcat(oPath, filesep, 'ERP_C_non_EBL_all.mat')), 'ERP_C_non_All')


    % Import TF analysis and ERP data by groups
    for k = 1:length(gr) % Loop through all groups
        % Initialize matrix for this groupa

        % Power_E_gr(:,:,:) = nan(num_chan_tf, num_frex, num_time);
        % Power_C_gr(:,:,:) = nan(num_chan_tf, num_frex, num_time);
        % ITPC_E_gr(:,:,:) =  nan(num_chan_tf, num_frex, num_time);
        % ITPC_C_gr(:,:,:) =  nan(num_chan_tf, num_frex, num_time);
        ERP_E_gr(:,:) = nan(num_chan_erp, num_time);
        ERP_C_gr(:,:) = nan(num_chan_erp, num_time);
        ERP_C_non_gr(:,:) = nan(num_chan_erp, num_time);

        for m = 1:num_chan_erp % Loop through all channels
            % Initialize matrix for this channel's data

            % Power_E_temp_gr = nan(length(Participants), num_frex, num_time);
            % Power_C_temp_gr = nan(length(Participants), num_frex, num_time);
            % ITPC_E_temp_gr = nan(length(Participants), num_frex, num_time);
            % ITPC_C_temp_gr = nan(length(Participants), num_frex, num_time);
            ERP_E_temp_gr = nan(length(Participants), num_time);
            ERP_C_temp_gr = nan(length(Participants), num_time);
            ERP_C_non_temp_gr = nan(length(Participants), num_time);

            for j = 1:length(Participants) % Loop through all participants for each channel
                load(char(strcat([Participants(j)], '_data_EBL.mat')))
                if strcmp(data.group,gr{k})
                    % Use data with more 8 or more usable errors
                    if data.UsableError > 7
                        % Power_E_temp_gr(j,:,:) = data.power.error(m,:,:);
                        % Power_C_temp_gr(j,:,:) = data.power.correct(m,:,:);
                        % ITPC_E_temp_gr(j,:,:) = data.ITPC.error(m,:,:);
                        % ITPC_C_temp_gr(j,:,:) = data.ITPC.correct(m,:,:);

                        % if m <= num_chan_erp % ERP channels have 2 less channels due to referencing to mastoid-equivalents
                        ERP_E_temp_gr(j,:) = data.erp.error(m,:);
                        ERP_C_temp_gr(j,:) = data.erp.correct(m,:);
                        ERP_C_non_temp_gr(j,:) = data.erp.correct_non(m,:);
                        % end % End ERP
                    end % End usable data extraction
                end % End participant loop for this particular channel
            end
            % Calculate mean across participants for this channel
            % Power_E_gr(m,:,:) = squeeze(mean(Power_E_temp_gr(:,:,:), 1, 'omitnan'));
            % Power_C_gr(m,:,:) = squeeze(mean(Power_C_temp_gr(:,:,:), 1, 'omitnan'));
            % ITPC_E_gr(m,:,:) =  squeeze(mean(ITPC_E_temp_gr(:,:,:), 1, 'omitnan'));
            % ITPC_C_gr(m,:,:) =  squeeze(mean(ITPC_C_temp_gr(:,:,:), 1, 'omitnan'));
            % if m <= num_chan_erp
            ERP_E_gr(m,:) = squeeze(mean(ERP_E_temp_gr(:,:), 1, 'omitnan'));
            ERP_C_gr(m,:) = squeeze(mean(ERP_C_temp_gr(:,:), 1, 'omitnan'));
            ERP_C_non_gr(m,:) = squeeze(mean(ERP_C_non_temp_gr(:,:), 1, 'omitnan'));
            % end
        end

        % Save data with average for this group
        % save(char(strcat(oPath, filesep, 'Power_E_', gr{k}, '.mat')), 'Power_E_gr')
        % save(char(strcat(oPath, filesep, 'Power_C_', gr{k}, '.mat')), 'Power_C_gr')
        % save(char(strcat(oPath, filesep, 'ITPC_E_', gr{k}, '.mat')), 'ITPC_E_gr')
        % save(char(strcat(oPath, filesep, 'ITPC_C_', gr{k}, '.mat')), 'ITPC_C_gr')
        save(char(strcat(oPath, filesep, 'ERP_E_EBL_', gr{k}, '.mat')), 'ERP_E_gr')
        save(char(strcat(oPath, filesep, 'ERP_C_EBL_', gr{k}, '.mat')), 'ERP_C_gr')
        save(char(strcat(oPath, filesep, 'ERP_C_non_EBL_', gr{k}, '.mat')), 'ERP_C_non_gr')
    end
end


%% Step 4 equivalent

for fol = 1:length(folders)
    % Power_E_All = [];
    % Power_E_BD = [];
    % Power_E_HC = [];
    % Power_E_SZ = [];
    % Power_C_All = [];
    % Power_C_BD = [];
    % Power_C_HC = [];
    % Power_C_SZ = [];
    % Power_diff = [];
    % ITPC_E_All = [];
    % ITPC_E_BD = [];
    % ITPC_E_HC = [];
    % ITPC_E_SZ = [];
    % ITPC_C_All = [];
    % ITPC_C_BD = [];
    % ITPC_C_HC = [];
    % ITPC_C_SZ = [];
    ERP_E_All = [];
    ERP_E_BD = [];
    ERP_E_HC = [];
    ERP_E_SZ = [];
    ERP_C_All = [];
    ERP_C_BD = [];
    ERP_C_HC = [];
    ERP_C_SZ = [];
    ERP_C_non_All = [];
    ERP_C_non_BD = [];
    ERP_C_non_HC = [];
    ERP_C_non_SZ = [];
    ERP_diff = [];
    ERP_diff_non = [];
    Analysis_Data = [];
    Data_Export = [];
    temp = [];


    % Define paths and create folders
    TaskPath = [ParentPath filesep folders{fol} filesep 'Output_Auto_updated'];
    iPath = [TaskPath filesep 'Response_PowerITPC_EBL'];
    oPath = [TaskPath filesep 'AnalysisData_EBL'];

    % Load files with "_data.mat" suffix in the
    cd(iPath);
    files = dir('*_data_EBL.mat');
    filenames = {files.name};
    Participants = regexprep(filenames, '_data_EBL.mat', '');

    load(char(strcat([Participants(1)], '_data_EBL.mat')))

    %% Time-Frequnecy plots
    % Load data with average from all participants to identify maximums
    % This part is not necessary, but just in case something happened and
    % restarting from the middle is helpful
    % load(char(strcat(oPath, filesep, 'Power_E_all.mat')))
    % load(char(strcat(oPath, filesep, 'Power_C_all.mat')))
    % load(char(strcat(oPath, filesep, 'ITPC_E_all.mat')))
    % load(char(strcat(oPath, filesep, 'ITPC_C_all.mat')))
    load(char(strcat(oPath, filesep, 'ERP_E_EBL_all.mat')))
    load(char(strcat(oPath, filesep, 'ERP_C_EBL_all.mat')))
    load(char(strcat(oPath, filesep, 'ERP_C_non_EBL_all.mat')))

    for g = 1:length(gr)
        % Power_E_gr = [];
        % Power_C_gr = [];
        % ITPC_E_gr = [];
        % ITPC_C_gr = [];
        ERP_E_gr = [];
        ERP_C_gr = [];
        ERP_C_non_gr = [];

        % load(char(strcat(oPath, filesep, 'Power_E_', gr{g}, '.mat')))
        % load(char(strcat(oPath, filesep, 'Power_C_', gr{g}, '.mat')))
        % load(char(strcat(oPath, filesep, 'ITPC_E_', gr{g}, '.mat')))
        % load(char(strcat(oPath, filesep, 'ITPC_C_', gr{g}, '.mat')))
        load(char(strcat(oPath, filesep, 'ERP_E_EBL_', gr{g}, '.mat')))
        load(char(strcat(oPath, filesep, 'ERP_C_EBL_', gr{g}, '.mat')))
        load(char(strcat(oPath, filesep, 'ERP_C_non_EBL_', gr{g}, '.mat')))

        if g == 1
            % Rename HC data
            % Power_E_HC = Power_E_gr;
            % Power_C_HC = Power_C_gr;
            % ITPC_E_HC = ITPC_E_gr;
            % ITPC_C_HC = ITPC_C_gr;
            ERP_E_HC = ERP_E_gr;
            ERP_C_HC = ERP_C_gr;
            ERP_C_non_HC = ERP_C_non_gr;

        elseif g == 2
            % Power_E_SZ = Power_E_gr;
            % Power_C_SZ = Power_C_gr;
            % ITPC_E_SZ = ITPC_E_gr;
            % ITPC_C_SZ = ITPC_C_gr;
            ERP_E_SZ = ERP_E_gr;
            ERP_C_SZ = ERP_C_gr;
            ERP_C_non_SZ = ERP_C_non_gr;

        elseif g == 3
            % Power_E_BD = Power_E_gr;
            % Power_C_BD = Power_C_gr;
            % ITPC_E_BD = ITPC_E_gr;
            % ITPC_C_BD = ITPC_C_gr;
            ERP_E_BD = ERP_E_gr;
            ERP_C_BD = ERP_C_gr;
            ERP_C_non_BD = ERP_C_non_gr;
        end
    end


    if fol == 1
        % Power_E_All_Arr = Power_E_All;
        % Power_C_All_Arr = Power_C_All;
        % ITPC_E_All_Arr = ITPC_E_All;
        % ITPC_C_All_Arr = ITPC_C_All;
        ERP_E_All_Arr = ERP_E_All;
        ERP_C_All_Arr = ERP_C_All;
        ERP_C_non_All_Arr = ERP_C_non_All;

        % Power_E_HC_Arr = Power_E_HC;
        % Power_C_HC_Arr = Power_C_HC;
        % ITPC_E_HC_Arr = ITPC_E_HC;
        % ITPC_C_HC_Arr = ITPC_C_HC;
        ERP_E_HC_Arr = ERP_E_HC;
        ERP_C_HC_Arr = ERP_C_HC;
        ERP_C_non_HC_Arr = ERP_C_non_HC;

        % Power_E_SZ_Arr = Power_E_SZ;
        % Power_C_SZ_Arr = Power_C_SZ;
        % ITPC_E_SZ_Arr = ITPC_E_SZ;
        % ITPC_C_SZ_Arr = ITPC_C_SZ;
        ERP_E_SZ_Arr = ERP_E_SZ;
        ERP_C_SZ_Arr = ERP_C_SZ;
        ERP_C_non_SZ_Arr = ERP_C_non_SZ;

        % Power_E_BD_Arr = Power_E_BD;
        % Power_C_BD_Arr = Power_C_BD;
        % ITPC_E_BD_Arr = ITPC_E_BD;
        % ITPC_C_BD_Arr = ITPC_C_BD;
        ERP_E_BD_Arr = ERP_E_BD;
        ERP_C_BD_Arr = ERP_C_BD;
        ERP_C_non_BD_Arr = ERP_C_non_BD;

    elseif fol == 2
        % Power_E_All_Neg = Power_E_All;
        % Power_C_All_Neg = Power_C_All;
        % ITPC_E_All_Neg = ITPC_E_All;
        % ITPC_C_All_Neg = ITPC_C_All;
        ERP_E_All_Neg = ERP_E_All;
        ERP_C_All_Neg = ERP_C_All;
        ERP_C_non_All_Neg = ERP_C_non_All;

        % Power_E_HC_Neg = Power_E_HC;
        % Power_C_HC_Neg = Power_C_HC;
        % ITPC_E_HC_Neg = ITPC_E_HC;
        % ITPC_C_HC_Neg = ITPC_C_HC;
        ERP_E_HC_Neg = ERP_E_HC;
        ERP_C_HC_Neg = ERP_C_HC;
        ERP_C_non_HC_Neg = ERP_C_non_HC;

        % Power_E_SZ_Neg = Power_E_SZ;
        % Power_C_SZ_Neg = Power_C_SZ;
        % ITPC_E_SZ_Neg = ITPC_E_SZ;
        % ITPC_C_SZ_Neg = ITPC_C_SZ;
        ERP_E_SZ_Neg = ERP_E_SZ;
        ERP_C_SZ_Neg = ERP_C_SZ;
        ERP_C_non_SZ_Neg = ERP_C_non_SZ;

        % Power_E_BD_Neg = Power_E_BD;
        % Power_C_BD_Neg = Power_C_BD;
        % ITPC_E_BD_Neg = ITPC_E_BD;
        % ITPC_C_BD_Neg = ITPC_C_BD;
        ERP_E_BD_Neg = ERP_E_BD;
        ERP_C_BD_Neg = ERP_C_BD;
        ERP_C_non_BD_Neg = ERP_C_non_BD;

    elseif fol == 3
        % Power_E_All_Pos = Power_E_All;
        % Power_C_All_Pos = Power_C_All;
        % ITPC_E_All_Pos = ITPC_E_All;
        % ITPC_C_All_Pos = ITPC_C_All;
        ERP_E_All_Pos = ERP_E_All;
        ERP_C_All_Pos = ERP_C_All;
        ERP_C_non_All_Pos = ERP_C_non_All;

        % Power_E_HC_Pos = Power_E_HC;
        % Power_C_HC_Pos = Power_C_HC;
        % ITPC_E_HC_Pos = ITPC_E_HC;
        % ITPC_C_HC_Pos = ITPC_C_HC;
        ERP_E_HC_Pos = ERP_E_HC;
        ERP_C_HC_Pos = ERP_C_HC;
        ERP_C_non_HC_Pos = ERP_C_non_HC;

        % Power_E_SZ_Pos = Power_E_SZ;
        % Power_C_SZ_Pos = Power_C_SZ;
        % ITPC_E_SZ_Pos = ITPC_E_SZ;
        % ITPC_C_SZ_Pos = ITPC_C_SZ;
        ERP_E_SZ_Pos = ERP_E_SZ;
        ERP_C_SZ_Pos = ERP_C_SZ;
        ERP_C_non_SZ_Pos = ERP_C_non_SZ;

        % Power_E_BD_Pos = Power_E_BD;
        % Power_C_BD_Pos = Power_C_BD;
        % ITPC_E_BD_Pos = ITPC_E_BD;
        % ITPC_C_BD_Pos = ITPC_C_BD;
        ERP_E_BD_Pos = ERP_E_BD;
        ERP_C_BD_Pos = ERP_C_BD;
        ERP_C_non_BD_Pos = ERP_C_non_BD;

    end
end


%% Find maximum and minim difference area (only within contrast, ignore
% participants)

% Start with power difference calculations
% Power_diff = Power_E_All - Power_C_All;
ERP_E_All_Combined = cat(3, ERP_E_All_Arr, ERP_E_All_Neg, ERP_E_All_Pos);
ERP_E_All_Combined_Average = squeeze(mean(ERP_E_All_Combined(:,:,:), 3, 'omitnan'));
ERP_C_All_Combined = cat(3, ERP_C_All_Arr, ERP_C_All_Neg, ERP_C_All_Pos);
ERP_C_All_Combined_Average = squeeze(mean(ERP_C_All_Combined(:,:,:), 3, 'omitnan'));
ERP_Diff_All_Combined = ERP_E_All_Combined_Average-ERP_C_All_Combined_Average;

ERP_C_non_All_Combined = cat(3, ERP_C_non_All_Arr, ERP_C_non_All_Neg, ERP_C_non_All_Pos);
ERP_C_non_All_Combined_Average = squeeze(mean(ERP_C_non_All_Combined(:,:,:), 3, 'omitnan'));
ERP_Diff_non_All_Combined = ERP_E_All_Combined_Average-ERP_C_non_All_Combined_Average;


% THETA at Cz (Global max)
% [max_tf_num,max_tf_idx] = max(Power_diff(:));
% [max_tf_l,max_tf_f,max_tf_t] = ind2sub(size(Power_diff),max_tf_idx);
%
% chan_tf_max = data.chan_tf(max_tf_l).labels;
% chan_n_tf_theta = find(ismember({data.chan_tf.labels}, chan_tf_max));
%
% freq_tf_max_cen = data.freq(max_tf_f);
% freq_tf_max_lower = data.freq(max_tf_f) - 2;
% freq_tf_max_width = 4;
% freq_tf_max_upper = freq_tf_max_lower + freq_tf_max_width;
%
% time_tf_max_cen = data.time(max_tf_t);
% time_tf_max_lower = data.time(max_tf_t) - 100;
% time_tf_max_width = 200;
% time_tf_max_upper = time_tf_max_lower + time_tf_max_width;
%
% % Find theta indices to extract data later
% freq_tf_max_lower_idx = find(data.freq <= freq_tf_max_lower, 1, 'last');
% freq_tf_max_upper_idx = find(data.freq >= freq_tf_max_upper, 1, 'first');
% time_tf_max_lower_idx = find(data.time <= time_tf_max_lower, 1, 'last');
% time_tf_max_upper_idx = find(data.time >= time_tf_max_upper, 1, 'first');

% Find maximum time of maximum ERN and save indices (matched)
% Note: ERN is "negative" so finding maximum ERN requires finding
% minimum value
chan_n_ern = find(ismember({data.chan_erp.labels}, 'Cz'));
[max_ERN_num, max_ERN_idx] = min(ERP_Diff_All_Combined(chan_n_ern, :));

ERN_cen = data.time(max_ERN_idx);
ERN_lower = data.time(max_ERN_idx) - 50;
ERN_width = 100;
ERN_upper = ERN_lower + ERN_width;

% Find indices for later
ERN_lower_idx = find(data.time <= ERN_lower, 1, 'last');
ERN_upper_idx = find(data.time >= ERN_upper, 1, 'first');

% Find maximum time of maximum ERN and save indices (non-matched)
[max_ERN_non_num, max_ERN_non_idx] = min(ERP_Diff_non_All_Combined(chan_n_ern, :));

ERN_non_cen = data.time(max_ERN_non_idx);
ERN_non_lower = data.time(max_ERN_non_idx) - 50;
ERN_non_width = 100;
ERN_non_upper = ERN_non_lower + ERN_non_width;

% Find indices for later
ERN_non_lower_idx = find(data.time <= ERN_non_lower, 1, 'last');
ERN_non_upper_idx = find(data.time >= ERN_non_upper, 1, 'first');

cPath = [ParentPath filesep 'Flanker_Outputs_updated' filesep 'Output_Auto_Combined'];

if ~isdir ([cPath])
    mkdir ([cPath]);
end


% Initialize data-driven max/min export table
max_export = table(ERN_cen, ERN_non_cen);
max_export.Properties.VariableNames = {'ERN_cen', 'ERN_non_cen'};
writetable(max_export, [cPath filesep 'Automatic_Parameters_EBL_' datestr(now,'yyyy-mm-dd'),'.csv']);

for fol = 1:length(folders)
    % Define paths and create folders
    TaskPath = [ParentPath filesep folders{fol} filesep 'Output_Auto_updated'];
    iPath = [TaskPath filesep 'Response_PowerITPC_EBL'];

    % Load files with "_data.mat" suffix in the
    cd(iPath);
    files = dir('*_data_EBL.mat');
    filenames = {files.name};
    Participants = regexprep(filenames, '_data_EBL.mat', '');
    %% Exporting all data at defined area from each task
    Data_Export = zeros(length(Participants), 10);

    for i = 1:length(Participants)
        cd(iPath);
        load([Participants{i}, '_data_EBL.mat']);

        % %% Extracting only Theta
        % %Save theta error power
        % temp = data.power.error(chan_n_tf_theta, freq_tf_max_lower_idx:freq_tf_max_upper_idx, time_tf_max_lower_idx:time_tf_max_upper_idx);
        % Data_Export(i, 1) = mean(temp(:));
        % %Save theta correct power
        % temp = data.power.correct(chan_n_tf_theta, freq_tf_max_lower_idx:freq_tf_max_upper_idx, time_tf_max_lower_idx:time_tf_max_upper_idx);
        % Data_Export(i, 2) = mean(temp(:));
        % %Save theta error itpc
        % temp = data.ITPC.error(chan_n_tf_theta, freq_tf_max_lower_idx:freq_tf_max_upper_idx, time_tf_max_lower_idx:time_tf_max_upper_idx);
        % Data_Export(i, 3) = mean(temp(:));
        % %Save theta correct itpc
        % temp = data.ITPC.correct(chan_n_tf_theta ,freq_tf_max_lower_idx:freq_tf_max_upper_idx, time_tf_max_lower_idx:time_tf_max_upper_idx);
        % Data_Export(i, 4) = mean(temp(:));

        %Save ERN
        temp = data.erp.error(chan_n_ern, ERN_lower_idx:ERN_upper_idx);
        Data_Export(i, 1) = mean(temp(:));
        %Save CRN trials matched to ERN
        temp = data.erp.correct(chan_n_ern, ERN_lower_idx:ERN_upper_idx);
        Data_Export(i, 2) = mean(temp(:));
        %Save ERN based on all CRN trials
        temp = data.erp.error(chan_n_ern, ERN_non_lower_idx:ERN_non_upper_idx);
        Data_Export(i, 3) = mean(temp(:));
        %Save all CRN trials
        temp = data.erp.correct_non(chan_n_ern, ERN_non_lower_idx:ERN_non_upper_idx);
        Data_Export(i, 4) = mean(temp(:));

        %Add other info
        Data_Export(i, 5) = data.UsableError;
        Data_Export(i, 6) = data.UsableCorrect;

        Data_Export(i,7) = data.AllError;
        Data_Export(i,8) = data.AllCorrect;

        Data_Export(i,9) = data.accuracy;

        if data.group == "HC"
            Data_Export(i,10) = 1;
        elseif data.group == "SZ"
            Data_Export(i,10) = 2;
        elseif data.group == "BD"
            Data_Export(i,10) = 3;
        else
            Data_Export(i,10) = NA;
        end
    end

    % Initialize data to export
    Analysis_Data = nan(length(Participants), (size(Data_Export, 2)+1));

    % Combine data
    Analysis_Data = [Participants', num2cell(Data_Export)];
    Analysis_Data = array2table(Analysis_Data);
    Analysis_Data.Properties.VariableNames = ["ID", "ERN","CRN", "ERN_non", "CRN_non",...
        "Error_n", "Correct_n", "Error_all", "Correct_all", "Accuracy", "Group"];

    % Below "writematrix" is preferred, but for version/software reasons you
    % cannot, use "save" option (2nd line, currently muted)
    writetable(Analysis_Data, char([cPath filesep 'EEG_' folders{fol} '_Data_Resp_EBL_' datestr(now,'yyyy-mm-dd') '.csv']));
    save(char(strcat([cPath filesep 'EEG_' folders{fol} '_Data_Resp_EBL_' datestr(now,'yyyy-mm-dd') '.mat'])), 'Analysis_Data');
end
params = [];

% Save parameters used in figures for each task
% params.chan_tf_max = chan_tf_max;
% params.chan_n_tf_theta = chan_n_tf_theta ;
% params.freq_tf_max_lower = freq_tf_max_lower;
% params.freq_tf_max_width = freq_tf_max_width;
% params.time_tf_max_lower = time_tf_max_lower ;
% params.time_tf_max_width = time_tf_max_width ;
% params.freq_tf_max_lower_idx = freq_tf_max_lower_idx ;
% params.freq_tf_max_upper_idx = freq_tf_max_upper_idx ;
% params.time_tf_max_lower_idx = time_tf_max_lower_idx ;
% params.time_tf_max_upper_idx = time_tf_max_upper_idx ;
params.chan_n_ern = chan_n_ern;
params.ERN_lower = ERN_lower ;
params.ERN_upper = ERN_upper ;
params.ERN_lower_idx = ERN_lower_idx ;
params.ERN_upper_idx = ERN_upper_idx ;
params.ERN_non_lower = ERN_non_lower;
params.ERN_non_upper = ERN_non_upper ;
params.ERN_non_lower_idx = ERN_non_lower_idx;
params.ERN_non_upper_idx = ERN_non_upper_idx ;

save(char(strcat(oPath, filesep, 'parameters_EBL.mat')), 'params');




fPath = [ParentPath filesep 'Flanker_Outputs_updated' filesep 'Figures_EBL'];

if ~isdir ([fPath])
    mkdir ([fPath]);
end

% ERN - matched
set(0,'defaultAxesFontSize',14, 'defaultTextFontSize',14)

figure('Position', [100 100 1100 600])
tiles = tiledlayout(2, 3,'TileSpacing','compact','Padding','compact');

% All Error and Correct bys
nexttile(1)
plot(data.time,squeeze(ERP_E_All_Arr(chan_n_ern,:)), data.time, squeeze(ERP_C_All_Arr(chan_n_ern,:)))
set(gca, 'xlim', [-400, 600], 'ylim', [-12, 12])
patch([params.ERN_lower params.ERN_lower params.ERN_upper params.ERN_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
ylabel({'Error and Correct'; 'Amplitude (\muV)'}, 'FontWeight', 'bold')
title('Arrow')

nexttile(2)
plot(data.time,squeeze(ERP_E_All_Neg(chan_n_ern,:)), data.time, squeeze(ERP_C_All_Neg(chan_n_ern,:)))
set(gca, 'xlim', [-400, 600], 'ylim', [-12, 12])
patch([params.ERN_lower params.ERN_lower params.ERN_upper params.ERN_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
title('Negative')

nexttile(3)
plot(data.time,squeeze(ERP_E_All_Pos(chan_n_ern,:)), data.time, squeeze(ERP_C_All_Pos(chan_n_ern,:)))
set(gca, 'xlim', [-400, 600], 'ylim', [-12, 12])
patch([params.ERN_lower params.ERN_lower params.ERN_upper params.ERN_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
legend('Error', 'Correct', 'Location', 'southeast')
title('Positive')

% All Differences bys
nexttile(4)
plot(data.time,squeeze(ERP_E_All_Arr(chan_n_ern,:)-ERP_C_All_Arr(chan_n_ern,:)), "k")
set(gca, 'xlim', [-400, 600], 'ylim', [-12, 12])
patch([params.ERN_lower params.ERN_lower params.ERN_upper params.ERN_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
ylabel({'Difference'; 'Amplitude (\muV)'}, 'FontWeight', 'bold')
xlabel('Time (ms)', 'FontWeight', 'bold')

nexttile(5)
plot(data.time,squeeze(ERP_E_All_Neg(chan_n_ern,:)-ERP_C_All_Neg(chan_n_ern,:)), "k")
set(gca, 'xlim', [-400, 600], 'ylim', [-12, 12])
patch([params.ERN_lower params.ERN_lower params.ERN_upper params.ERN_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
xlabel('Time (ms)', 'FontWeight', 'bold')

nexttile(6)
plot(data.time,squeeze(ERP_E_All_Pos(chan_n_ern,:)-ERP_C_All_Pos(chan_n_ern,:)), "k")
set(gca, 'xlim', [-400, 600], 'ylim', [-12, 12])
patch([params.ERN_lower params.ERN_lower params.ERN_upper params.ERN_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
legend('Difference')
xlabel('Time (ms)', 'FontWeight', 'bold')

title(tiles, 'ERN, CRN, and \DeltaERN', 'FontSize', 16, 'FontWeight', 'Bold')

saveas(gcf, [fPath filesep 'Resp_ERN_All_EBL_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_ERN_All_EBL_' datestr(now,'yyyy-mm-dd') '.fig'])



% ERN – Non-Matched
figure('Position', [100 100 1100 600])
tiles = tiledlayout(2, 3,'TileSpacing','compact','Padding','compact');

% All Error and Correct bys
nexttile(1)
plot(data.time,squeeze(ERP_E_All_Arr(chan_n_ern,:)), data.time, squeeze(ERP_C_non_All_Arr(chan_n_ern,:)))
set(gca, 'xlim', [-400, 600], 'ylim', [-12, 12])
patch([params.ERN_non_lower params.ERN_non_lower params.ERN_non_upper params.ERN_non_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
ylabel({'Error and Correct'; 'Amplitude (\muV)'}, 'FontWeight', 'bold')
title('Arrow')

nexttile(2)
plot(data.time,squeeze(ERP_E_All_Neg(chan_n_ern,:)), data.time, squeeze(ERP_C_non_All_Neg(chan_n_ern,:)))
set(gca, 'xlim', [-400, 600], 'ylim', [-12, 12])
patch([params.ERN_non_lower params.ERN_non_lower params.ERN_non_upper params.ERN_non_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
title('Negative')

nexttile(3)
plot(data.time,squeeze(ERP_E_All_Pos(chan_n_ern,:)), data.time, squeeze(ERP_C_non_All_Pos(chan_n_ern,:)))
set(gca, 'xlim', [-400, 600], 'ylim', [-12, 12])
patch([params.ERN_non_lower params.ERN_non_lower params.ERN_non_upper params.ERN_non_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
legend('Error', 'Correct (All Trials)', 'Location', 'southeast')
title('Positive')

% All Differences bys
nexttile(4)
plot(data.time,squeeze(ERP_E_All_Arr(chan_n_ern,:)-ERP_C_non_All_Arr(chan_n_ern,:)), "k")
set(gca, 'xlim', [-400, 600], 'ylim', [-12, 12])
patch([params.ERN_non_lower params.ERN_non_lower params.ERN_non_upper params.ERN_non_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
ylabel({'Difference'; 'Amplitude (\muV)'}, 'FontWeight', 'bold')

nexttile(5)
plot(data.time,squeeze(ERP_E_All_Neg(chan_n_ern,:)-ERP_C_non_All_Neg(chan_n_ern,:)), "k")
set(gca, 'xlim', [-400, 600], 'ylim', [-12, 12])
patch([params.ERN_non_lower params.ERN_non_lower params.ERN_non_upper params.ERN_non_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])

nexttile(6)
plot(data.time,squeeze(ERP_E_All_Pos(chan_n_ern,:)-ERP_C_non_All_Pos(chan_n_ern,:)), "k")
set(gca, 'xlim', [-400, 600], 'ylim', [-12, 12])
patch([params.ERN_non_lower params.ERN_non_lower params.ERN_non_upper params.ERN_non_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
legend('Difference')

title(tiles, 'ERN, CRN (All Trials), and \DeltaERN', 'FontSize', 16, 'FontWeight', 'Bold')

saveas(gcf, [fPath filesep 'Resp_ERN_non_All_EBL_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_ERN_non_All_EBL_' datestr(now,'yyyy-mm-dd') '.fig'])


% ERN topography
figure('Position', [100 100 1100 850])

limit_Arr = max(abs(mean(ERP_C_All_Arr(:, params.ERN_lower_idx: params.ERN_upper_idx), 2, 'omitnan')));
limit_Neg = max(abs(mean(ERP_C_All_Neg(:, params.ERN_lower_idx: params.ERN_upper_idx), 2, 'omitnan')));
limit_Pos = max(abs(mean(ERP_C_All_Pos(:, params.ERN_lower_idx: params.ERN_upper_idx), 2, 'omitnan')));
limit = max([limit_Arr, limit_Neg, limit_Pos]);
topo_limits = [-limit limit];

% All Error bys
subplot(3,3,1)
topoplot((mean(ERP_E_All_Arr(:, params.ERN_lower_idx: params.ERN_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
text(-.8, 0, 'Error', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 14)
colormap(jet)
title('Arrow', 'FontSize', 14)

subplot(3,3,2)
topoplot((mean(ERP_E_All_Neg(:, params.ERN_lower_idx: params.ERN_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)
title('Negative', 'FontSize', 14)

subplot(3,3,3)
topoplot((mean(ERP_E_All_Pos(:, params.ERN_lower_idx: params.ERN_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)
title('Positive', 'FontSize', 14)

c = colorbar('eastoutside');
pos = get(gca, 'Position');
c.Position = [pos(1) + pos(3) + .01, pos(2), .02, pos(4)];

% All Correct bys
subplot(3,3,4)
topoplot((mean(ERP_C_All_Arr(:, params.ERN_lower_idx: params.ERN_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
text(-.8, 0, 'Correct', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 14)
colormap(jet)

subplot(3,3,5)
topoplot((mean(ERP_C_All_Neg(:, params.ERN_lower_idx: params.ERN_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

subplot(3,3,6)
topoplot((mean(ERP_C_All_Pos(:, params.ERN_lower_idx: params.ERN_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

c = colorbar('eastoutside');
pos = get(gca, 'Position');
c.Position = [pos(1) + pos(3) + .01, pos(2), .02, pos(4)];

% All differences bys
subplot(3,3,7)
topoplot((mean(ERP_E_All_Arr(:, params.ERN_lower_idx: params.ERN_upper_idx), 2, 'omitnan') - mean(ERP_C_All_Arr(:, params.ERN_lower_idx: params.ERN_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
text(-.8, 0, 'Difference', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 14)
colormap(jet)

subplot(3,3,8)
topoplot((mean(ERP_E_All_Neg(:, params.ERN_lower_idx: params.ERN_upper_idx), 2, 'omitnan') - mean(ERP_C_All_Neg(:, params.ERN_lower_idx: params.ERN_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

subplot(3,3,9)
topoplot((mean(ERP_E_All_Pos(:, params.ERN_lower_idx: params.ERN_upper_idx), 2, 'omitnan') - mean(ERP_C_All_Pos(:, params.ERN_lower_idx: params.ERN_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

c = colorbar('eastoutside');
pos = get(gca, 'Position');
c.Position = [pos(1) + pos(3) + .01, pos(2), .02, pos(4)];

annotation('textbox', [.52, 1, 0, 0], 'String', 'ERN, CRN, and \DeltaERN', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'Bold', 'FitBoxToText', 'on', 'EdgeColor', 'none')

saveas(gcf, [fPath filesep 'Resp_ERN_Topography_All_EBL_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_ERN_Topography_All_EBL_' datestr(now,'yyyy-mm-dd') '.fig'])



% ERN Non-Matched topography
figure('Position', [100 100 1100 850])

limit_Arr = max(abs(mean(ERP_C_non_All_Arr(:, params.ERN_non_lower_idx: params.ERN_non_upper_idx), 2, 'omitnan')));
limit_Neg = max(abs(mean(ERP_C_non_All_Neg(:, params.ERN_non_lower_idx: params.ERN_non_upper_idx), 2, 'omitnan')));
limit_Pos = max(abs(mean(ERP_C_non_All_Pos(:, params.ERN_non_lower_idx: params.ERN_non_upper_idx), 2, 'omitnan')));
limit = max([limit_Arr, limit_Neg, limit_Pos]);
topo_limits = [-limit limit];

% All Error bys
subplot(3,3,1)
topoplot((mean(ERP_E_All_Arr(:, params.ERN_non_lower_idx: params.ERN_non_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
text(-.8, 0, 'Error', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 14)
colormap(jet)
title('Arrow', 'FontSize', 14)

subplot(3,3,2)
topoplot((mean(ERP_E_All_Neg(:, params.ERN_non_lower_idx: params.ERN_non_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)
title('Negative', 'FontSize', 14)

subplot(3,3,3)
topoplot((mean(ERP_E_All_Pos(:, params.ERN_non_lower_idx: params.ERN_non_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)
title('Positive', 'FontSize', 14)

c = colorbar('eastoutside');
pos = get(gca, 'Position');
c.Position = [pos(1) + pos(3) + .01, pos(2), .02, pos(4)];

% All Correct bys
subplot(3,3,4)
topoplot((mean(ERP_C_non_All_Arr(:, params.ERN_non_lower_idx: params.ERN_non_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
text(-.8, 0, 'Correct', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 14)
colormap(jet)

subplot(3,3,5)
topoplot((mean(ERP_C_non_All_Neg(:, params.ERN_non_lower_idx: params.ERN_non_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

subplot(3,3,6)
topoplot((mean(ERP_C_non_All_Pos(:, params.ERN_non_lower_idx: params.ERN_non_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

c = colorbar('eastoutside');
pos = get(gca, 'Position');
c.Position = [pos(1) + pos(3) + .01, pos(2), .02, pos(4)];

% All Differences bys
subplot(3,3,7)
topoplot((mean(ERP_E_All_Arr(:, params.ERN_non_lower_idx: params.ERN_non_upper_idx), 2, 'omitnan') - mean(ERP_C_non_All_Arr(:, params.ERN_non_lower_idx: params.ERN_non_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
text(-.8, 0, 'Difference', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 14)
colormap(jet)

subplot(3,3,8)
topoplot((mean(ERP_E_All_Neg(:, params.ERN_non_lower_idx: params.ERN_non_upper_idx), 2, 'omitnan') - mean(ERP_C_non_All_Neg(:, params.ERN_non_lower_idx: params.ERN_non_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

subplot(3,3,9)
topoplot((mean(ERP_E_All_Pos(:, params.ERN_non_lower_idx: params.ERN_non_upper_idx), 2, 'omitnan') - mean(ERP_C_non_All_Pos(:, params.ERN_non_lower_idx: params.ERN_non_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

c = colorbar('eastoutside');
pos = get(gca, 'Position');
c.Position = [pos(1) + pos(3) + .01, pos(2), .02, pos(4)];

annotation('textbox', [.52, 1, 0, 0], 'String', 'ERN, CRN (All Trials), and \DeltaERN', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'Bold', 'FitBoxToText', 'on', 'EdgeColor', 'none')

saveas(gcf, [fPath filesep 'Resp_Theta_ERP_non_Topography_All_EBL_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_Theta_ERP_non_Topography_All_EBL_' datestr(now,'yyyy-mm-dd') '.fig'])



fprintf('Step 5 alternative baseline calculations are done \n');