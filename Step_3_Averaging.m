%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Response Monitoring Theta-Band Activities across Emotional Contexts in 
% Schizophrenia- and Bipolar-Spectrum Disorders
% Suzuki, Menkes, et al.
%
% Script to conduct compile all data extracted in Step 2
%
% Completed using MATLAB 2024b & EEGLAB 2024.0, on Windows 11 Enterprise
%
% Author: Takakuni Suzuki
% First drafted June 2023
% Last updated January 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
clear;
clc;


%%%%%%%%% Setting up folders %%%%%%%%%
ParentPath = ''; % Specifcy data folder
folders = {'FARR' 'FNEG' 'FPOS'};
gr = {'HC','SZ','BD'};


for fol = 1:length(folders)
    TaskPath = [ParentPath filesep folders{fol} filesep 'Output_Auto_updated'];

    %% Define paths and create folders
    iPath = [TaskPath filesep 'Response_PowerITPC'];
    oPath = [TaskPath filesep 'AnalysisData'];

    if ~isdir ([oPath])
        mkdir ([oPath]);
    end

    % Load files with "_data.mat" suffix
    cd(iPath);
    files = dir('*_data.mat');
    filenames = {files.name};
    Participants = regexprep(filenames, '_data.mat', '');

    % Import first data to extract parameters (e.g., channel locations)
    load(char(strcat([Participants(1)], '_data.mat')))
    num_chan_tf = length(data.chan_tf);
    num_chan_erp = length(data.chan_erp);
    num_frex = length(data.freq);
    num_time = length(data.time);

    % Initiate TF and ERP data matrices for all participants
    Power_E_All(:,:,:) = nan(num_chan_tf, num_frex, num_time);
    Power_C_All(:,:,:) = nan(num_chan_tf, num_frex, num_time);
    ITPC_E_All(:,:,:) =  nan(num_chan_tf, num_frex, num_time);
    ITPC_C_All(:,:,:) =  nan(num_chan_tf, num_frex, num_time);
    ERP_E_All(:,:) = nan(num_chan_erp, num_time);
    ERP_C_All(:,:) = nan(num_chan_erp, num_time);
    ERP_C_non_All(:,:) = nan(num_chan_erp, num_time);

    %% Import TF and ERP data from all participants
    % Going through channels to avoid computer memory issues
    for l = 1:num_chan_tf % loop through all channels
        % Initialize matrix for this channel's data
        Power_E_All_temp = nan(length(Participants), num_frex, num_time);
        Power_C_All_temp = nan(length(Participants), num_frex, num_time);
        ITPC_E_All_temp = nan(length(Participants), num_frex, num_time);
        ITPC_C_All_temp = nan(length(Participants), num_frex, num_time);
        ERP_E_All_temp = nan(length(Participants), num_time);
        ERP_C_All_temp = nan(length(Participants), num_time);
        ERP_C_All_non_temp = nan(length(Participants),  num_time);

        for i = 1:length(Participants) % Loop through all participants for each channel
            load(char(strcat([Participants(i)], '_data.mat')))

            % Use data with 8 or more usable errors
            if data.UsableError > 7
                Power_E_All_temp(i,:,:) = data.power.error(l,:,:);
                Power_C_All_temp(i,:,:) = data.power.correct(l,:,:);
                ITPC_E_All_temp(i,:,:) = data.ITPC.error(l,:,:);
                ITPC_C_All_temp(i,:,:) = data.ITPC.correct(l,:,:);

                if l <= num_chan_erp % ERP channels have 2 less channels due to referencing to mastoid-equivalents
                    ERP_E_All_temp(i,:) = data.erp.error(l,:);
                    ERP_C_All_temp(i,:) = data.erp.correct(l,:);
                    ERP_C_All_non_temp(i,:) = data.erp.correct_non(l,:);
                end % End ERP
            end % End usable data extraction
        end % End participant loop for this particular channel

        % Calculate mean across participants for this channel
        Power_E_All(l,:,:) = squeeze(mean(Power_E_All_temp(:,:,:), 1, 'omitnan'));
        Power_C_All(l,:,:) = squeeze(mean(Power_C_All_temp(:,:,:), 1, 'omitnan'));
        ITPC_E_All(l,:,:) =  squeeze(mean(ITPC_E_All_temp(:,:,:), 1, 'omitnan'));
        ITPC_C_All(l,:,:) =  squeeze(mean(ITPC_C_All_temp(:,:,:), 1, 'omitnan'));

        if l <= num_chan_erp
            ERP_E_All(l,:) = squeeze(mean(ERP_E_All_temp(:,:), 1, 'omitnan'));
            ERP_C_All(l,:)  = squeeze(mean(ERP_C_All_temp(:,:), 1, 'omitnan'));
            ERP_C_non_All(l,:)  = squeeze(mean(ERP_C_All_non_temp(:,:), 1, 'omitnan'));
        end
    end

    % Save data with average from all participants
    save(char(strcat(oPath, filesep, 'Power_E_all.mat')), 'Power_E_All')
    save(char(strcat(oPath, filesep, 'Power_C_all.mat')), 'Power_C_All')
    save(char(strcat(oPath, filesep, 'ITPC_E_all.mat')), 'ITPC_E_All')
    save(char(strcat(oPath, filesep, 'ITPC_C_all.mat')), 'ITPC_C_All')
    save(char(strcat(oPath, filesep, 'ERP_E_all.mat')), 'ERP_E_All')
    save(char(strcat(oPath, filesep, 'ERP_C_all.mat')), 'ERP_C_All')
    save(char(strcat(oPath, filesep, 'ERP_C_non_all.mat')), 'ERP_C_non_All')


    % Import TF analysis and ERP data by groups
    for k = 1:length(gr) % Loop through all groups
        % Initialize matrix for this groupa

        Power_E_gr(:,:,:) = nan(num_chan_tf, num_frex, num_time);
        Power_C_gr(:,:,:) = nan(num_chan_tf, num_frex, num_time);
        ITPC_E_gr(:,:,:) =  nan(num_chan_tf, num_frex, num_time);
        ITPC_C_gr(:,:,:) =  nan(num_chan_tf, num_frex, num_time);
        ERP_E_gr(:,:) = nan(num_chan_erp, num_time);
        ERP_C_gr(:,:) = nan(num_chan_erp, num_time);
        ERP_C_non_gr(:,:) = nan(num_chan_erp, num_time);

        for m = 1:num_chan_tf % Loop through all channels
            % Initialize matrix for this channel's data

            Power_E_temp_gr = nan(length(Participants), num_frex, num_time);
            Power_C_temp_gr = nan(length(Participants), num_frex, num_time);
            ITPC_E_temp_gr = nan(length(Participants), num_frex, num_time);
            ITPC_C_temp_gr = nan(length(Participants), num_frex, num_time);
            ERP_E_temp_gr = nan(length(Participants), num_time);
            ERP_C_temp_gr = nan(length(Participants), num_time);
            ERP_C_non_temp_gr = nan(length(Participants), num_time);

            for j = 1:length(Participants) % Loop through all participants for each channel
                load(char(strcat([Participants(j)], '_data.mat')))
                if strcmp(data.group,gr{k})
                    % Use data with 8 or more usable errors
                    if data.UsableError > 7 && data.accuracy > .65
                        Power_E_temp_gr(j,:,:) = data.power.error(m,:,:);
                        Power_C_temp_gr(j,:,:) = data.power.correct(m,:,:);
                        ITPC_E_temp_gr(j,:,:) = data.ITPC.error(m,:,:);
                        ITPC_C_temp_gr(j,:,:) = data.ITPC.correct(m,:,:);

                        if m <= num_chan_erp % ERP channels have 2 less channels due to referencing to mastoid-equivalents
                            ERP_E_temp_gr(j,:) = data.erp.error(m,:);
                            ERP_C_temp_gr(j,:) = data.erp.correct(m,:);
                            ERP_C_non_temp_gr(j,:) = data.erp.correct_non(m,:);
                        end % End ERP
                    end % End usable data extraction
                end % End participant loop for this particular channel
            end
            % Calculate mean across participants for this channel
            Power_E_gr(m,:,:) = squeeze(mean(Power_E_temp_gr(:,:,:), 1, 'omitnan'));
            Power_C_gr(m,:,:) = squeeze(mean(Power_C_temp_gr(:,:,:), 1, 'omitnan'));
            ITPC_E_gr(m,:,:) =  squeeze(mean(ITPC_E_temp_gr(:,:,:), 1, 'omitnan'));
            ITPC_C_gr(m,:,:) =  squeeze(mean(ITPC_C_temp_gr(:,:,:), 1, 'omitnan'));
            if m <= num_chan_erp
                ERP_E_gr(m,:) = squeeze(mean(ERP_E_temp_gr(:,:), 1, 'omitnan'));
                ERP_C_gr(m,:) = squeeze(mean(ERP_C_temp_gr(:,:), 1, 'omitnan'));
                ERP_C_non_gr(m,:) = squeeze(mean(ERP_C_non_temp_gr(:,:), 1, 'omitnan'));
            end
        end

        % Save data with average for this group
        save(char(strcat(oPath, filesep, 'Power_E_', gr{k}, '.mat')), 'Power_E_gr')
        save(char(strcat(oPath, filesep, 'Power_C_', gr{k}, '.mat')), 'Power_C_gr')
        save(char(strcat(oPath, filesep, 'ITPC_E_', gr{k}, '.mat')), 'ITPC_E_gr')
        save(char(strcat(oPath, filesep, 'ITPC_C_', gr{k}, '.mat')), 'ITPC_C_gr')
        save(char(strcat(oPath, filesep, 'ERP_E_', gr{k}, '.mat')), 'ERP_E_gr')
        save(char(strcat(oPath, filesep, 'ERP_C_', gr{k}, '.mat')), 'ERP_C_gr')
        save(char(strcat(oPath, filesep, 'ERP_C_non_', gr{k}, '.mat')), 'ERP_C_non_gr')
    end
end
fprintf('Step 3 Average calculations are done \n');