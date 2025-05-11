%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Response Monitoring Theta-Band Activities across Emotional Contexts in 
% Schizophrenia- and Bipolar-Spectrum Disorders
% Suzuki, Menkes, et al.
%
% Script to extract time-frequency and ERP from data preprocessed in Step 1
%
% Codes based on scripts shared by Cohen, M. X. (2014). Analyzing neural 
% time series data: theory and practice. MIT press.
% Mostly from Ch. 13
%
% Completed using MATLAB 2024b & EEGLAB 2024.0, on Windows 11 Enterprise
%
% Author: Takakuni Suzuki
% First drafted June 2023
% Last updated May 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

%%%%%%%%% Setting up the preprocessing parameters, folders, etc. %%%%%%%%%
ParentPath = ''; % Specifcy data folder
folders = {'FARR' 'FNEG' 'FPOS'};

% Define frequency range to analyze
min_freq = 2; % We're intersted starting around 3 hz-ish. This is to have frequencies in power of 2
max_freq = 128;
num_frex = 43;

for fol = 1:length(folders)
    % Set paths
    TaskPath = [ParentPath filesep folders{fol}];
    iPath = [TaskPath filesep 'Preprocess_Auto' filesep '3_processed_response'];
    oParentPath = [TaskPath filesep 'Output_Auto_updated'];
    o1Path = [oParentPath filesep 'Response_Matched'];
    o2Path = [oParentPath filesep 'Response_PowerITPC'];

    % Set suffix
    isuffix = '_processed_response.set';
    osuffix = '_matched_response.set';

    % Enter response marker info
    Resp_Cor = {'S  3' 'S  4'};
    Resp_Err = {'S  9' 'S 10'};
    Resp_Events = [Resp_Cor, Resp_Err];
    Stim_Events = {'S  1' 'S  2'};

    window = [-3 (-3 + (500 * 5.5 - 1)/500)];

    % Create folders
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
        row = [];
        EEG = [];
        EEG_base = [];
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

        % Initialize time frequency data
        basepower = NaN(EEG.nbchan,num_frex);
        eegpower_error = NaN(EEG.nbchan,num_frex,EEG.pnts);
        eegpower_correct= NaN(EEG.nbchan,num_frex,EEG.pnts);
        eegitpc_error = NaN(EEG.nbchan,num_frex,EEG.pnts);
        eegitpc_correct = NaN(EEG.nbchan,num_frex,EEG.pnts);
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

            %% Removed
            % This was an error on the 1st submission and pre-print. 
            % This nullifies the trial-matching, but retaining it here for
            % transparency

            % % Re-read data to re-import stimuli and re-epoch
            % EEG = pop_loadset([Participants_n{i}, isuffix], iPath);
            % EEG = pop_epoch( EEG, Resp_Events, window);
            %
            % % To remvoe event codes from previous and after the responses
            % EEG = pop_selectevent( EEG, 'latency','-100 <= 1','type', {'S  3', 'S  4', 'S  9', 'S 10'},...
            % 'deleteevents','on','deleteepochs','off','invertepochs','off');
            %% End removed section 

            % Reject unneeded epochs
            EEG = pop_selectevent( EEG, 'epoch', idxEpo ,'deleteevents','off','deleteepochs','on','invertepochs','off');

            % Save the data that only keeps correct trials matched to error trials.
            pop_saveset(EEG, [Participants_n{i} osuffix], o1Path);

            %% Start time-frequency decomposition
            %% TF of baseline
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

            data.n_trials_all = length(EEG_base.epoch);

            % Define wavelet parameters
            time = -1:1/EEG_base.srate:1;
            frex = logspace(log10(min_freq),log10(max_freq),num_frex);
            s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);

            % Define convolution parameters
            n_wavelet_base            = length(time);
            n_data_base               = EEG_base.pnts*EEG_base.trials;
            n_convolution_base        = n_wavelet_base+n_data_base-1;
            n_conv_pow2_base          = pow2(nextpow2(n_convolution_base));
            half_of_wavelet_size_base = (n_wavelet_base-1)/2;

            % Loop through each electrode
            for k = 1:EEG_base.nbchan
                chan2use_base = EEG_base.chanlocs(k).labels;

                % Get FFT of data
                eegfft_base = fft(reshape(EEG_base.data(strcmpi(chan2use_base,{EEG_base.chanlocs.labels}),:,:),1,EEG_base.pnts*EEG_base.trials),n_conv_pow2_base);

                % Find locations of the baseline window rangea
                baseidx = dsearchn(EEG_base.times',[-300 -50]');

                % Loop through frequencies and compute synchronization
                for fi=1:num_frex
                    wavelet_base = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2_base );

                    % Convolution
                    eegconv_base = ifft(wavelet_base.*eegfft_base);
                    eegconv_base = eegconv_base(1:n_convolution_base);
                    eegconv_base = eegconv_base(half_of_wavelet_size_base+1:end-half_of_wavelet_size_base);

                    % Average power over trials
                    temppower_base = mean(abs(reshape(eegconv_base,EEG_base.pnts,EEG_base.trials)).^2,2);
                    basepower(k,fi,:) = mean(temppower_base(baseidx(1):baseidx(2)));
                end
            end

            %% TF of error trials
            EEG_error = pop_loadset([Participants_n{i} osuffix], o1Path);

            % Extract only error trials
            EEG_error = pop_selectevent( EEG_error, 'type',Resp_Err,'deleteevents','off','deleteepochs','on','invertepochs','off');
            EEG_error = eeg_checkset( EEG_error );

            data.n_trials_error = length(EEG_error.epoch);

            % Define wavelet parameters
            time = -1:1/EEG_error.srate:1;
            frex = logspace(log10(min_freq),log10(max_freq),num_frex);
            s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);

            % Define convolution parameters
            n_wavelet_error            = length(time);
            n_data_error               = EEG_error.pnts*EEG_error.trials;
            n_convolution_error        = n_wavelet_error+n_data_error-1;
            n_conv_pow2_error          = pow2(nextpow2(n_convolution_error));
            half_of_wavelet_size_error = (n_wavelet_error-1)/2;

            % Loop through each electrode
            for j = 1:EEG_error.nbchan %EEG.nbchan
                chan2use_error = EEG_error.chanlocs(j).labels;

                % Get FFT of data
                eegfft_error = fft(reshape(EEG_error.data(strcmpi(chan2use_error,{EEG_error.chanlocs.labels}),:,:),1,EEG_error.pnts*EEG_error.trials),n_conv_pow2_error);

                % Loop through frequencies and compute synchronization
                for fi=1:num_frex
                    wavelet_error = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2_error );

                    % convolution
                    eegconv_error = ifft(wavelet_error.*eegfft_error);
                    eegconv_error = eegconv_error(1:n_convolution_error);
                    eegconv_error = eegconv_error(half_of_wavelet_size_error+1:end-half_of_wavelet_size_error);

                    % Average power over trials and adjust for baseline
                    temppower_error = mean(abs(reshape(eegconv_error, EEG_error.pnts, EEG_error.trials)).^2,2);
                    eegpower_error(j,fi,:) = 10*log10(temppower_error./basepower(j,fi));
                    eegitpc_error(j,fi,:) = abs(mean(exp(1i*angle(reshape(eegconv_error,EEG_error.pnts,EEG_error.trials))),2));
                end
            end
            tf_chan = EEG_error.chanlocs;


            % For ERN
            EEG_error = pop_reref( EEG_error, [10 21] );% Mastoid (tp9/tp10) reference
            EEG_error = pop_rmbase( EEG_error, [-200 -50]); % Baseline
            for k = 1:EEG_error.nbchan %EEG.nbchan
                erp_error(k,:) = nanmean(EEG_error.data(k,:,:),3) ;
            end
            erp_chan = EEG_error.chanlocs;

            %% TF for Correct trials
            EEG_correct = pop_loadset([Participants_n{i} osuffix], o1Path);

            % Define wavelet parameters
            time = -1:1/EEG_correct.srate:1;
            frex = logspace(log10(min_freq),log10(max_freq),num_frex);
            s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);

            % Extract only correct trials
            EEG_correct = pop_selectevent( EEG_correct, 'type',Resp_Cor,'deleteevents','off','deleteepochs','on','invertepochs','off');
            EEG_correct = eeg_checkset( EEG_correct );

            data.n_trials_correct= length(EEG_correct.epoch);

            % Define convolution parameters
            n_wavelet_correct            = length(time);
            n_data_correct               = EEG_correct.pnts*EEG_correct.trials;
            n_convolution_correct        = n_wavelet_correct+n_data_correct-1;
            n_conv_pow2_correct          = pow2(nextpow2(n_convolution_correct));
            half_of_wavelet_size_correct = (n_wavelet_correct-1)/2;

            % Loop through each electrode
            for l = 1:EEG_correct.nbchan %EEG.nbchan
                chan2use_correct = EEG_correct.chanlocs(l).labels;

                % Get FFT of data
                eegfft_correct = fft(reshape(EEG_correct.data(strcmpi(chan2use_correct,{EEG_correct.chanlocs.labels}),:,:),1,EEG_correct.pnts*EEG_correct.trials),n_conv_pow2_correct);

                % Loop through frequencies and compute synchronization
                for fi=1:num_frex
                    wavelet_correct = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2_correct );

                    % Convolution
                    eegconv_correct = ifft(wavelet_correct.*eegfft_correct);
                    eegconv_correct = eegconv_correct(1:n_convolution_correct);
                    eegconv_correct = eegconv_correct(half_of_wavelet_size_correct+1:end-half_of_wavelet_size_correct);

                    % Average power over trials and adjust for baseline
                    temppower_correct = mean(abs(reshape(eegconv_correct,EEG_correct.pnts,EEG_correct.trials)).^2,2);
                    eegpower_correct(l,fi,:) = 10*log10(temppower_correct./basepower(l,fi));
                    eegitpc_correct(l,fi,:) = abs(mean(exp(1i*angle(reshape(eegconv_correct,EEG_correct.pnts,EEG_correct.trials))),2));
                end
            end

            % For matched CRN
            EEG_correct = pop_reref( EEG_correct, [10 21] );% Mastoid (tp9/tp10) reference
            EEG_correct = pop_rmbase( EEG_correct, [-200 -50]); %Baseline
            for m = 1:EEG_correct.nbchan %EEG.nbchan
                erp_correct(m,:) = nanmean(EEG_correct.data(m,:,:),3) ;
            end

            % For non-matched CRN
            EEG_correct_non = pop_loadset([Participants_n{i}, isuffix], iPath);
            EEG_correct_non = pop_epoch( EEG_correct_non, Resp_Cor, window);
            EEG_correct_non = pop_reref( EEG_correct_non, [10 21] );% Mastoid (tp9/tp10) reference
            EEG_correct_non = pop_rmbase( EEG_correct_non, [-200 -50]); % Baseline
            for n = 1:EEG_correct_non.nbchan
                erp_correct_non(n,:) = nanmean(EEG_correct_non.data(n,:,:),3) ;
            end

            data.power.error = eegpower_error;
            data.power.correct = eegpower_correct;
            data.ITPC.error = eegitpc_error;
            data.ITPC.correct = eegitpc_correct;
            data.erp.error = erp_error;
            data.erp.correct = erp_correct;
            data.erp.correct_non = erp_correct_non;
            data.freq = frex;
            data.time = EEG.times;
            data.chan_tf = tf_chan;
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

            save(char(strcat(o2Path, filesep, Participants_n(i), '_data', '.mat')), 'data');
        end
    end
end
fprintf('Step 2 TF and ERP extraction are done \n');
