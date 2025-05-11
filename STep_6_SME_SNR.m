%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Response Monitoring Theta-Band Activities across Emotional Contexts in 
% Schizophrenia- and Bipolar-Spectrum Disorders
% Suzuki, Menkes, et al.
%
% Script to calculate standardized measurement error and
% signal-to-noise ratio 
%
% Methods described in
% Luck, Stewart, Simmons, and Rhemtulla (2021). Standardized measurement
% error: A universal metric of data quality for averaaged event-related
% potentials. Psychophysiology, 58:e13793.
% https://doi.org/10.1111/psyp.13793
%
% Completed using MATLAB 2024b & EEGLAB 2024.0, on Windows 11 Enterprise
%
% Author: Takakuni Suzuki
% First drafted December, 2024
% Last updated May 2025
%
% This was a step suggested by anonymous reviewers. The second script could
% be adjusted to make this process more efficient (i.e., saving trial-level
% data at intermediate steps), but the scripts used for the manuscript is
% uploaded.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

%%%%%%%%% Setting up the preprocessing parameters, folders, etc. %%%%%%%%%
ParentPath = ''; % Specifcy data folder
folders = {'FARR' 'FNEG' 'FPOS'};
cPath = [ParentPath filesep 'Flanker_Outputs_updated\Output_Auto_Combined'];
parameters = readtable([cPath filesep 'Automatic_Parameters_2025-05-11.csv']);

%% Define frequency range to analyze
min_freq = 2; % We're intersted starting around 3 hz-ish. This is to have frequencies in power of 2 
max_freq = 128;
num_frex = 43;

for fol = 1:length(folders)
    SME_data = [];

    % Define paths and create folders
    TaskPath = [ParentPath filesep folders{fol}];
    iPath = [TaskPath filesep 'Preprocess_Auto' filesep '3_processed_response'];
    oParentPath = [TaskPath filesep 'Output_Auto_updated'];
    o1Path = [oParentPath filesep 'Response_Matched'];
    % o2Path = [oParentPath filesep 'Response_PowerITPC'];

    % Set suffix
    isuffix = '_processed_response.set';
    osuffix = '_matched_response.set';

    % Enter response marker info
    Resp_Cor = {'S  3' 'S  4'};
    Resp_Err = {'S  9' 'S 10'};
    Resp_Events = [Resp_Cor, Resp_Err];
    Stim_Events = {'S  1' 'S  2'};

    window = [-3 (-3 + (500 * 5.5 - 1)/500)];

    %%%%% Start analysis
    %% Read in participant files
    cd(iPath);
    Participants=dir(['*' isuffix]);
    Participants={Participants.name};

    % Remove extra parts of the file names
    Participants_n = string(nan(length(Participants), 1));
    for j = 1:length(Participants)
        name = char(Participants(j));
        Participants_n(j) = string(name(1:13));
    end

    SME_data = nan(length(Participants),9);

    %% Start with matching number of correct trails to number of error trials
    for i = 1:length(Participants_n)
        % Initialize data matrix
        data = [];
        row = [];
        EEG = [];
        EEG_base = [];
        EEG_correct= [];
        EEG_correct_non = [];
        EEG_error = [];

        eegpower_error_trials = [];
        eegpower_correct_trials = [];

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
        basepower = NaN(1,num_frex);
        eegpower_error = NaN(num_frex,EEG.pnts);
        eegpower_correct= NaN(num_frex,EEG.pnts);
        eegitpc_error = NaN(num_frex,EEG.pnts);
        eegitpc_correct = NaN(num_frex,EEG.pnts);
        erp_error = NaN(1,EEG.pnts);
        erp_correct = NaN(1,EEG.pnts);
        erp_correct_non = NaN(1,EEG.pnts);


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

            %% Start Time Frequency Analyses
            % Exract Baseline
            EEG_base = pop_loadset([Participants_n{i}, isuffix], iPath);

            idx_base = setdiff([EEG_base.event.urevent], idxEpo_base);
            for h = 1:length(EEG_base.event)
                if  ismember(EEG_base.event(h).type, Stim_Events) && ismember(EEG_base.event(h).urevent, idx_base)
                    EEG_base.event(h).type = '98';
                end
            end

            % Extract epochs based on stimulus presentation
            % Extracting 300ms extra on the front end for baseline window
            EEG_base = pop_epoch( EEG_base, Stim_Events, [-1.8 1.51]);
            EEG_base = eeg_checkset( EEG_base );

            % Define wavelet parameters
            time = -1:1/EEG_base.srate:1;
            frex = logspace(log10(min_freq),log10(max_freq),num_frex);
            s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);

            % Definte convolution parameters
            n_wavelet_base            = length(time);
            n_data_base               = EEG_base.pnts*EEG_base.trials;
            n_convolution_base        = n_wavelet_base+n_data_base-1;
            n_conv_pow2_base          = pow2(nextpow2(n_convolution_base));
            half_of_wavelet_size_base = (n_wavelet_base-1)/2;

            % Select channel
            chan_n_tf_theta = find(ismember({EEG_base.chanlocs.labels}, parameters.chan_tf_max));

            % get FFT of data
            eegfft_base = fft(reshape(EEG_base.data(chan_n_tf_theta,:,:),1,EEG_base.pnts*EEG_base.trials),n_conv_pow2_base);

            % Find locations of the baseline window range
            baseidx = dsearchn(EEG_base.times',[-300 -50]');

            % loop through frequencies and compute synchronization
            for fi=1:num_frex
                wavelet_base = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2_base );

                % convolution
                eegconv_base = ifft(wavelet_base.*eegfft_base);
                eegconv_base = eegconv_base(1:n_convolution_base);
                eegconv_base = eegconv_base(half_of_wavelet_size_base+1:end-half_of_wavelet_size_base);

                % Average power over trials
                temppower_base = mean(abs(reshape(eegconv_base,EEG_base.pnts,EEG_base.trials)).^2,2);
                basepower(fi) = mean(temppower_base(baseidx(1):baseidx(2)));
            end


            %% TF for error trials
            EEG_error = pop_loadset([Participants_n{i} osuffix], o1Path);

            % extract only error trials
            EEG_error = pop_selectevent( EEG_error, 'type',Resp_Err,'deleteevents','off','deleteepochs','on','invertepochs','off');
            EEG_error = eeg_checkset( EEG_error );

            % define wavelet parameters
            time = -1:1/EEG_error.srate:1;
            frex = logspace(log10(min_freq),log10(max_freq),num_frex);
            s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);

            % definte convolution parameters
            n_wavelet_error            = length(time);
            n_data_error               = EEG_error.pnts*EEG_error.trials;
            n_convolution_error        = n_wavelet_error+n_data_error-1;
            n_conv_pow2_error          = pow2(nextpow2(n_convolution_error));
            half_of_wavelet_size_error = (n_wavelet_error-1)/2;

            % Define extraction window
            freq_tf_max_lower = parameters.freq_tf_max_cen - 2;
            freq_tf_max_width = 4;
            freq_tf_max_upper = freq_tf_max_lower + freq_tf_max_width;

            time_tf_max_lower = parameters.time_tf_max_cen - 100;
            time_tf_max_width = 200;
            time_tf_max_upper = time_tf_max_lower + time_tf_max_width;

            time_erp_min_lower = parameters.ERN_cen - 50;
            time_erp_min_width = 100;
            time_erp_min_upper = time_erp_min_lower + time_erp_min_width;

            % Find indices to extract data
            freq_tf_max_lower_idx = find(frex <= freq_tf_max_lower, 1, 'last');
            freq_tf_max_upper_idx = find(frex >= freq_tf_max_upper, 1, 'first');
            time_tf_max_lower_idx = find(EEG.times <= time_tf_max_lower, 1, 'last');
            time_tf_max_upper_idx = find(EEG.times >= time_tf_max_upper, 1, 'first');
            time_erp_min_lower_idx = find(EEG.times <= time_erp_min_lower, 1, 'last');
            time_erp_min_upper_idx = find(EEG.times >= time_erp_min_upper, 1, 'first');


            % get FFT of data
            eegfft_error = fft(reshape(EEG_error.data(chan_n_tf_theta,:,:),1,EEG_error.pnts*EEG_error.trials),n_conv_pow2_error);

            % loop through frequencies and compute synchronization
            for fi=1:num_frex
                wavelet_error = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2_error );

                % convolution
                eegconv_error = ifft(wavelet_error.*eegfft_error);
                eegconv_error = eegconv_error(1:n_convolution_error);
                eegconv_error = eegconv_error(half_of_wavelet_size_error+1:end-half_of_wavelet_size_error);

                % Adjust for baseline
                eegpower_error_trials_temp = 10*log10((abs(eegconv_error).^2)./basepower(fi));
                eegpower_error_trials(fi,:,:) = reshape(eegpower_error_trials_temp, EEG_error.pnts, EEG_error.trials);
            end

            sme_error_tf_temp = squeeze(mean(eegpower_error_trials(freq_tf_max_lower_idx:freq_tf_max_upper_idx, time_tf_max_lower_idx:time_tf_max_upper_idx,:),1:2));
            SME_data(i,1) = std(sme_error_tf_temp)/sqrt(length(sme_error_tf_temp));

            %For ERN
            EEG_error = pop_reref( EEG_error, [10 21] );% Mastoid (tp9/tp10) reference
            EEG_error = pop_rmbase( EEG_error, [-200 -50]); %Baseline
            chan_n_ern = find(ismember({EEG_error.chanlocs.labels}, parameters.chan_tf_max));

            erp_error_trials = squeeze(mean(EEG_error.data(chan_n_ern,time_erp_min_lower_idx:time_erp_min_upper_idx,:),1:2));
            SME_data(i,2) = std(erp_error_trials)/sqrt(length(erp_error_trials));
            SME_data(i,3) = length(sme_error_tf_temp);



            %% TF for Correct trials
            EEG_correct = pop_loadset([Participants_n{i} osuffix], o1Path);

            % define wavelet parameters
            time = -1:1/EEG_correct.srate:1;
            frex = logspace(log10(min_freq),log10(max_freq),num_frex);
            s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);

            % extract only correct trials
            EEG_correct = pop_selectevent( EEG_correct, 'type',Resp_Cor,'deleteevents','off','deleteepochs','on','invertepochs','off');
            EEG_correct = eeg_checkset( EEG_correct );

            % definte convolution parameters
            n_wavelet_correct            = length(time);
            n_data_correct               = EEG_correct.pnts*EEG_correct.trials;
            n_convolution_correct        = n_wavelet_correct+n_data_correct-1;
            n_conv_pow2_correct          = pow2(nextpow2(n_convolution_correct));
            half_of_wavelet_size_correct = (n_wavelet_correct-1)/2;

            % get FFT of data
            eegfft_correct = fft(reshape(EEG_correct.data(chan_n_tf_theta,:,:),1,EEG_correct.pnts*EEG_correct.trials),n_conv_pow2_correct);

            % loop through frequencies and compute synchronization
            for fi=1:num_frex
                wavelet_correct = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2_correct );

                % convolution
                eegconv_correct = ifft(wavelet_correct.*eegfft_correct);
                eegconv_correct = eegconv_correct(1:n_convolution_correct);
                eegconv_correct = eegconv_correct(half_of_wavelet_size_correct+1:end-half_of_wavelet_size_correct);

                % Adjust for baseline
                eegpower_correct_trials_temp = 10*log10((abs(eegconv_correct).^2)./basepower(fi));
                eegpower_correct_trials(fi,:,:) = reshape(eegpower_correct_trials_temp, EEG_correct.pnts, EEG_correct.trials);

            end

            sme_correct_tf_temp = squeeze(mean(eegpower_correct_trials(freq_tf_max_lower_idx:freq_tf_max_upper_idx, time_tf_max_lower_idx:time_tf_max_upper_idx,:),1:2));
            SME_data(i,4) = std(sme_correct_tf_temp)/sqrt(length(sme_correct_tf_temp));

            %For matched CRN
            EEG_correct = pop_reref( EEG_correct, [10 21] );% Mastoid (tp9/tp10) reference
            EEG_correct = pop_rmbase( EEG_correct, [-200 -50]); %Baseline
            erp_correct_trials = squeeze(mean(EEG_correct.data(chan_n_ern,time_erp_min_lower_idx:time_erp_min_upper_idx,:),1:2));
            SME_data(i,5) = std(erp_correct_trials)/sqrt(length(erp_correct_trials));
            SME_data(i,6) = length(sme_correct_tf_temp);


            %For non-matched CRN
            EEG_correct_non = pop_loadset([Participants_n{i}, isuffix], iPath);
            EEG_correct_non = pop_epoch( EEG_correct_non, Resp_Cor, window);
            EEG_correct_non = eeg_checkset( EEG_correct_non );
            EEG_correct_non = pop_reref( EEG_correct_non, [10 21] );% Mastoid (tp9/tp10) reference
            EEG_correct_non = pop_rmbase( EEG_correct_non, [-200 -50]); %Baseline

            erp_correct_non_trials = squeeze(mean(EEG_correct_non.data(chan_n_ern,time_erp_min_lower_idx:time_erp_min_upper_idx,:),1:2));
            SME_data(i,7) = std(erp_correct_non_trials)/sqrt(length(erp_correct_non_trials));
            SME_data(i,8) = length(erp_correct_non_trials);

        end

        %Attaching group info
        if strcmp(Participants_n{i}(1:10),strcat('TNA_', folders{fol}, '_1'))
            SME_data(i,9) = 1;
        else
            if strcmp(Participants_n{i}(1:10),strcat('TNA_', folders{fol}, '_2'))
                SME_data(i,9) = 2;
            else
                if strcmp(Participants_n{i}(1:10),strcat('TNA_', folders{fol}, '_3'))
                    SME_data(i,9) = 3;
                else
                    SME_data(i,9) = [];
                end
            end
        end
    end

    SME_export = [Participants_n, num2cell(SME_data)];
    SME_export = array2table(SME_export);
    SME_export.Properties.VariableNames = {'ID', 'SME_theta_power_error', 'SME_ERN', 'Error_n', 'SME_theta_power_correct', 'SME_CRN', 'Correct_n', 'SME_CRN_non', 'Correct_non_n', 'Group'};
    writetable(SME_export, [cPath filesep 'SME_' folders{fol} '_' datestr(now,'yyyy-mm-dd'),'.csv']);
end

fprintf('Step 6 SME calculations are done \n');