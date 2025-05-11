%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Response Monitoring Theta-Band Activities across Emotional Contexts in 
% Schizophrenia- and Bipolar-Spectrum Disorders
% Suzuki, Menkes, et al.
%
% Script to create figures and extract data for statistical analyses
% from files created in Steps 2 and 3
%
% Codes based on scripts shared by Cohen, M. X. (2014). Analyzing neural 
% time series data: theory and practice. MIT press.
% particularly tfviewerx.m
%
% Process of finding max and min automated here, but checked manually and
% visually first
%
% Completed using MATLAB 2024b & EEGLAB 2024.0, on Windows 11 Enterprise
%
% Author: Takakuni Suzuki
% First drafted June 2023
% Last updated January 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

%%%%%%%%% Setting up the preprocessing parameters, folders, etc. %%%%%%%%%
ParentPath = ''; % Specifcy data folder
pPath = [ParentPath filesep 'Flanker_Outputs_updated'];
cPath = [pPath filesep 'Output_Auto_Combined'];
fPath = [pPath filesep 'Figures'];

folders = {'FARR' 'FNEG' 'FPOS'};
gr = {'HC','SZ','BD'};

if ~isdir ([cPath])
    mkdir ([cPath]);
end

if ~isdir ([fPath])
    mkdir ([fPath]);
end

% Frequency range to plot
min_freq = 2;
max_freq = 64; % Analyzed up to 128

eeglab

for fol = 1:length(folders)
    Power_E_All = [];
    Power_E_BD = [];
    Power_E_HC = [];
    Power_E_SZ = [];
    Power_C_All = [];
    Power_C_BD = [];
    Power_C_HC = [];
    Power_C_SZ = [];
    Power_diff = [];
    ITPC_E_All = [];
    ITPC_E_BD = [];
    ITPC_E_HC = [];
    ITPC_E_SZ = [];
    ITPC_C_All = [];
    ITPC_C_BD = [];
    ITPC_C_HC = [];
    ITPC_C_SZ = [];
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
    params = [];
    temp = [];

    % Define paths and create folders
    TaskPath = [ParentPath filesep folders{fol} filesep 'Output_Auto_updated'];
    iPath = [TaskPath filesep 'Response_PowerITPC'];
    oPath = [TaskPath filesep 'AnalysisData'];

    % Load files with "_data.mat" suffix in the
    cd(iPath);
    files = dir('*_data.mat');
    filenames = {files.name};
    Participants = regexprep(filenames, '_data.mat', '');

    load(char(strcat([Participants(1)], '_data.mat')))

    %% Time-Frequnecy plots
    % Load data with average from all participants to identify maximums
    load(char(strcat(oPath, filesep, 'Power_E_all.mat')))
    load(char(strcat(oPath, filesep, 'Power_C_all.mat')))
    load(char(strcat(oPath, filesep, 'ITPC_E_all.mat')))
    load(char(strcat(oPath, filesep, 'ITPC_C_all.mat')))
    load(char(strcat(oPath, filesep, 'ERP_E_all.mat')))
    load(char(strcat(oPath, filesep, 'ERP_C_all.mat')))
    load(char(strcat(oPath, filesep, 'ERP_C_non_all.mat')))

    for g = 1:length(gr)
        Power_E_gr = [];
        Power_C_gr = [];
        ITPC_E_gr = [];
        ITPC_C_gr = [];
        ERP_E_gr = [];
        ERP_C_gr = [];
        ERP_C_non_gr = [];

        load(char(strcat(oPath, filesep, 'Power_E_', gr{g}, '.mat')))
        load(char(strcat(oPath, filesep, 'Power_C_', gr{g}, '.mat')))
        load(char(strcat(oPath, filesep, 'ITPC_E_', gr{g}, '.mat')))
        load(char(strcat(oPath, filesep, 'ITPC_C_', gr{g}, '.mat')))
        load(char(strcat(oPath, filesep, 'ERP_E_', gr{g}, '.mat')))
        load(char(strcat(oPath, filesep, 'ERP_C_', gr{g}, '.mat')))
        load(char(strcat(oPath, filesep, 'ERP_C_non_', gr{g}, '.mat')))

        if g == 1
            % Rename HC data
            Power_E_HC = Power_E_gr;
            Power_C_HC = Power_C_gr;
            ITPC_E_HC = ITPC_E_gr;
            ITPC_C_HC = ITPC_C_gr;
            ERP_E_HC = ERP_E_gr;
            ERP_C_HC = ERP_C_gr;
            ERP_C_non_HC = ERP_C_non_gr;

        elseif g == 2
            Power_E_SZ = Power_E_gr;
            Power_C_SZ = Power_C_gr;
            ITPC_E_SZ = ITPC_E_gr;
            ITPC_C_SZ = ITPC_C_gr;
            ERP_E_SZ = ERP_E_gr;
            ERP_C_SZ = ERP_C_gr;
            ERP_C_non_SZ = ERP_C_non_gr;

        elseif g == 3
            Power_E_BD = Power_E_gr;
            Power_C_BD = Power_C_gr;
            ITPC_E_BD = ITPC_E_gr;
            ITPC_C_BD = ITPC_C_gr;
            ERP_E_BD = ERP_E_gr;
            ERP_C_BD = ERP_C_gr;
            ERP_C_non_BD = ERP_C_non_gr;
        end
    end


    if fol == 1
        Power_E_All_Arr = Power_E_All;
        Power_C_All_Arr = Power_C_All;
        ITPC_E_All_Arr = ITPC_E_All;
        ITPC_C_All_Arr = ITPC_C_All;
        ERP_E_All_Arr = ERP_E_All;
        ERP_C_All_Arr = ERP_C_All;
        ERP_C_non_All_Arr = ERP_C_non_All;

        Power_E_HC_Arr = Power_E_HC;
        Power_C_HC_Arr = Power_C_HC;
        ITPC_E_HC_Arr = ITPC_E_HC;
        ITPC_C_HC_Arr = ITPC_C_HC;
        ERP_E_HC_Arr = ERP_E_HC;
        ERP_C_HC_Arr = ERP_C_HC;
        ERP_C_non_HC_Arr = ERP_C_non_HC;

        Power_E_SZ_Arr = Power_E_SZ;
        Power_C_SZ_Arr = Power_C_SZ;
        ITPC_E_SZ_Arr = ITPC_E_SZ;
        ITPC_C_SZ_Arr = ITPC_C_SZ;
        ERP_E_SZ_Arr = ERP_E_SZ;
        ERP_C_SZ_Arr = ERP_C_SZ;
        ERP_C_non_SZ_Arr = ERP_C_non_SZ;

        Power_E_BD_Arr = Power_E_BD;
        Power_C_BD_Arr = Power_C_BD;
        ITPC_E_BD_Arr = ITPC_E_BD;
        ITPC_C_BD_Arr = ITPC_C_BD;
        ERP_E_BD_Arr = ERP_E_BD;
        ERP_C_BD_Arr = ERP_C_BD;
        ERP_C_non_BD_Arr = ERP_C_non_BD;

        params_Arr = params;

    elseif fol == 2
        Power_E_All_Neg = Power_E_All;
        Power_C_All_Neg = Power_C_All;
        ITPC_E_All_Neg = ITPC_E_All;
        ITPC_C_All_Neg = ITPC_C_All;
        ERP_E_All_Neg = ERP_E_All;
        ERP_C_All_Neg = ERP_C_All;
        ERP_C_non_All_Neg = ERP_C_non_All;

        Power_E_HC_Neg = Power_E_HC;
        Power_C_HC_Neg = Power_C_HC;
        ITPC_E_HC_Neg = ITPC_E_HC;
        ITPC_C_HC_Neg = ITPC_C_HC;
        ERP_E_HC_Neg = ERP_E_HC;
        ERP_C_HC_Neg = ERP_C_HC;
        ERP_C_non_HC_Neg = ERP_C_non_HC;

        Power_E_SZ_Neg = Power_E_SZ;
        Power_C_SZ_Neg = Power_C_SZ;
        ITPC_E_SZ_Neg = ITPC_E_SZ;
        ITPC_C_SZ_Neg = ITPC_C_SZ;
        ERP_E_SZ_Neg = ERP_E_SZ;
        ERP_C_SZ_Neg = ERP_C_SZ;
        ERP_C_non_SZ_Neg = ERP_C_non_SZ;

        Power_E_BD_Neg = Power_E_BD;
        Power_C_BD_Neg = Power_C_BD;
        ITPC_E_BD_Neg = ITPC_E_BD;
        ITPC_C_BD_Neg = ITPC_C_BD;
        ERP_E_BD_Neg = ERP_E_BD;
        ERP_C_BD_Neg = ERP_C_BD;
        ERP_C_non_BD_Neg = ERP_C_non_BD;

        params_Neg = params;

    elseif fol == 3
        Power_E_All_Pos = Power_E_All;
        Power_C_All_Pos = Power_C_All;
        ITPC_E_All_Pos = ITPC_E_All;
        ITPC_C_All_Pos = ITPC_C_All;
        ERP_E_All_Pos = ERP_E_All;
        ERP_C_All_Pos = ERP_C_All;
        ERP_C_non_All_Pos = ERP_C_non_All;

        Power_E_HC_Pos = Power_E_HC;
        Power_C_HC_Pos = Power_C_HC;
        ITPC_E_HC_Pos = ITPC_E_HC;
        ITPC_C_HC_Pos = ITPC_C_HC;
        ERP_E_HC_Pos = ERP_E_HC;
        ERP_C_HC_Pos = ERP_C_HC;
        ERP_C_non_HC_Pos = ERP_C_non_HC;

        Power_E_SZ_Pos = Power_E_SZ;
        Power_C_SZ_Pos = Power_C_SZ;
        ITPC_E_SZ_Pos = ITPC_E_SZ;
        ITPC_C_SZ_Pos = ITPC_C_SZ;
        ERP_E_SZ_Pos = ERP_E_SZ;
        ERP_C_SZ_Pos = ERP_C_SZ;
        ERP_C_non_SZ_Pos = ERP_C_non_SZ;

        Power_E_BD_Pos = Power_E_BD;
        Power_C_BD_Pos = Power_C_BD;
        ITPC_E_BD_Pos = ITPC_E_BD;
        ITPC_C_BD_Pos = ITPC_C_BD;
        ERP_E_BD_Pos = ERP_E_BD;
        ERP_C_BD_Pos = ERP_C_BD;
        ERP_C_non_BD_Pos = ERP_C_non_BD;

        params_Pos = params;
    end
end

Power_E_All_Combined = cat(4, Power_E_All_Arr, Power_E_All_Neg, Power_E_All_Pos);
Power_E_All_Combined_Average = squeeze(mean(Power_E_All_Combined(:,:,:,:), 4, 'omitnan'));
Power_C_All_Combined = cat(4, Power_C_All_Arr, Power_C_All_Neg, Power_C_All_Pos);
Power_C_All_Combined_Average = squeeze(mean(Power_C_All_Combined(:,:,:,:), 4, 'omitnan'));
Power_Diff_All_Combined = Power_E_All_Combined_Average-Power_C_All_Combined_Average;

% THETA at Cz (Global max)
[max_tf_num,max_tf_idx] = max(Power_Diff_All_Combined(:));
[max_tf_l,max_tf_f,max_tf_t] = ind2sub(size(Power_Diff_All_Combined),max_tf_idx);

chan_tf_max = data.chan_tf(max_tf_l).labels;
chan_n_tf_theta = find(ismember({data.chan_tf.labels}, chan_tf_max));

freq_tf_max_cen = data.freq(max_tf_f);
freq_tf_max_lower = data.freq(max_tf_f) - 2;
freq_tf_max_width = 4;
freq_tf_max_upper = freq_tf_max_lower + freq_tf_max_width;

time_tf_max_cen = data.time(max_tf_t);
time_tf_max_lower = data.time(max_tf_t) - 100;
time_tf_max_width = 200;
time_tf_max_upper = time_tf_max_lower + time_tf_max_width;

% Find theta indices to extract data later
freq_tf_max_lower_idx = find(data.freq <= freq_tf_max_lower, 1, 'last');
freq_tf_max_upper_idx = find(data.freq >= freq_tf_max_upper, 1, 'first');
time_tf_max_lower_idx = find(data.time <= time_tf_max_lower, 1, 'last');
time_tf_max_upper_idx = find(data.time >= time_tf_max_upper, 1, 'first');



ITPC_E_All_Combined = cat(4, ITPC_E_All_Arr, ITPC_E_All_Neg, ITPC_E_All_Pos);
ITPC_E_All_Combined_Average = squeeze(mean(ITPC_E_All_Combined(:,:,:,:), 4, 'omitnan'));
ITPC_C_All_Combined = cat(4, ITPC_C_All_Arr, ITPC_C_All_Neg, ITPC_C_All_Pos);
ITPC_C_All_Combined_Average = squeeze(mean(ITPC_C_All_Combined(:,:,:,:), 4, 'omitnan'));
ITPC_Diff_All_Combined = ITPC_E_All_Combined_Average-ITPC_C_All_Combined_Average;

ERP_E_All_Combined = cat(3, ERP_E_All_Arr, ERP_E_All_Neg, ERP_E_All_Pos);
ERP_E_All_Combined_Average = squeeze(mean(ERP_E_All_Combined(:,:,:), 3, 'omitnan'));
ERP_C_All_Combined = cat(3, ERP_C_All_Arr, ERP_C_All_Neg, ERP_C_All_Pos);
ERP_C_All_Combined_Average = squeeze(mean(ERP_C_All_Combined(:,:,:), 3, 'omitnan'));
ERP_Diff_All_Combined = ERP_E_All_Combined_Average-ERP_C_All_Combined_Average;

ERP_C_non_All_Combined = cat(3, ERP_C_non_All_Arr, ERP_C_non_All_Neg, ERP_C_non_All_Pos);
ERP_C_non_All_Combined_Average = squeeze(mean(ERP_C_non_All_Combined(:,:,:), 3, 'omitnan'));
ERP_Diff_non_All_Combined = ERP_E_All_Combined_Average-ERP_C_non_All_Combined_Average;

% Find maximum time of maximum ERN and save indices (matched)
% Note: ERN is "negative" so finding maximum ERN requires finding
% minimum value
chan_n_ern = find(ismember({data.chan_erp.labels}, chan_tf_max));
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

% Initialize data-driven max/min export table
max_export = table(convertCharsToStrings(chan_tf_max), freq_tf_max_cen, time_tf_max_cen, ERN_cen, ERN_non_cen);
max_export.Properties.VariableNames = {'chan_tf_max', 'freq_tf_max_cen', 'time_tf_max_cen', 'ERN_cen', 'ERN_non_cen'};
writetable(max_export, [cPath filesep 'Automatic_Parameters_' datestr(now,'yyyy-mm-dd'),'.csv']);

for fol = 1:length(folders)
    % Define paths and create folders
    TaskPath = [ParentPath filesep folders{fol} filesep 'Output_Auto_updated'];
    iPath = [TaskPath filesep 'Response_PowerITPC'];

    % Load files with "_data.mat" suffix in the
    cd(iPath);
    files = dir('*_data.mat');
    filenames = {files.name};
    Participants = regexprep(filenames, '_data.mat', '');

    %% Exporting all data at defined area from each task
    Data_Export = zeros(length(Participants), 17);

    for i = 1:length(Participants)
        cd(iPath);
        load([Participants{i}, '_data.mat']);

        % Extracting only Theta
        %Save theta error power
        temp = data.power.error(chan_n_tf_theta, freq_tf_max_lower_idx:freq_tf_max_upper_idx, time_tf_max_lower_idx:time_tf_max_upper_idx);
        Data_Export(i, 1) = mean(temp(:));
        %Save theta correct power
        temp = data.power.correct(chan_n_tf_theta, freq_tf_max_lower_idx:freq_tf_max_upper_idx, time_tf_max_lower_idx:time_tf_max_upper_idx);
        Data_Export(i, 2) = mean(temp(:));
        %Save theta error itpc
        temp = data.ITPC.error(chan_n_tf_theta, freq_tf_max_lower_idx:freq_tf_max_upper_idx, time_tf_max_lower_idx:time_tf_max_upper_idx);
        Data_Export(i, 3) = mean(temp(:));
        %Save theta correct itpc
        temp = data.ITPC.correct(chan_n_tf_theta ,freq_tf_max_lower_idx:freq_tf_max_upper_idx, time_tf_max_lower_idx:time_tf_max_upper_idx);
        Data_Export(i, 4) = mean(temp(:));

        %Save ERN
        temp = data.erp.error(chan_n_ern, ERN_lower_idx:ERN_upper_idx);
        Data_Export(i, 5) = mean(temp(:));
        %Save CRN trials matched to ERN
        temp = data.erp.correct(chan_n_ern, ERN_lower_idx:ERN_upper_idx);
        Data_Export(i, 6) = mean(temp(:));
        %Save ERN based on all CRN trials
        temp = data.erp.error(chan_n_ern, ERN_non_lower_idx:ERN_non_upper_idx);
        Data_Export(i, 7) = mean(temp(:));
        %Save all CRN trials
        temp = data.erp.correct_non(chan_n_ern, ERN_non_lower_idx:ERN_non_upper_idx);
        Data_Export(i, 8) = mean(temp(:));

        %Add other info
        Data_Export(i, 9) = data.UsableError;
        Data_Export(i, 10) = data.UsableCorrect;
        Data_Export(i,11) = data.AllError;
        Data_Export(i,12) = data.AllCorrect;
        Data_Export(i,13) = data.accuracy;

        if data.group == "HC"
            Data_Export(i,14) = 1;
        elseif data.group == "SZ"
            Data_Export(i,14) = 2;
        elseif data.group == "BD"
            Data_Export(i,14) = 3;
        else
            Data_Export(i,14) = NA;
        end
        
        Data_Export(i,15) = data.n_trials_all;
        Data_Export(i,16) = data.n_trials_error;
        Data_Export(i,17) = data.n_trials_correct;

    end

    % Initialize data to export
    Analysis_Data = nan(length(Participants), (size(Data_Export, 2)+1));

    % Combine data
    Analysis_Data = [Participants', num2cell(Data_Export)];
    Analysis_Data = array2table(Analysis_Data);
    Analysis_Data.Properties.VariableNames = ["ID", "Power_Theta_E", "Power_Theta_C", "ITPC_Theta_E", "ITPC_Theta_C", "ERN","CRN", "ERN_non", "CRN_non",...
        "Error_n", "Correct_n", "Error_all", "Correct_all", "Accuracy", "Group", "Matched_all_n", "Matched_error_n", "Matched_correct_n"];

    writetable(Analysis_Data, char([cPath filesep 'EEG_' folders{fol} '_Data_Resp_Combined_' datestr(now,'yyyy-mm-dd') '.csv']));
    save(char(strcat([cPath filesep 'EEG_' folders{fol} '_Data_Resp_Combined_' datestr(now,'yyyy-mm-dd') '.mat'])), 'Analysis_Data');
end 





%% Create figures
set(0,'defaultAxesFontSize',14, 'defaultTextFontSize',14)

% Theta
% All participants
% Power plots



% Max identification plot
figure('Position', [100 100 650 850])
tiles = tiledlayout(3, 1,'TileSpacing','compact','Padding','compact');

% All Error 
nexttile(1)
contourf(data.time,data.freq,squeeze(Power_E_All_Combined_Average(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Error'; 'Frequency (Hz)'}, 'FontWeight', 'bold')

% All Correct 
nexttile(2)
contourf(data.time,data.freq,squeeze(Power_C_All_Combined_Average(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Correct'; 'Frequency (Hz)'}, 'FontWeight', 'bold')


% Difference
nexttile(3)
contourf(data.time,data.freq,squeeze(Power_E_All_Combined_Average(chan_n_tf_theta,:,:) - Power_C_All_Combined_Average(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Difference'; 'Frequency (Hz)'}, 'FontWeight', 'bold')

title(tiles, 'Power (dB)', 'FontSize', 16, 'FontWeight', 'Bold')

saveas(gcf, [fPath filesep 'Resp_Theta_Power_All_Combined_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_Theta_Power_All_Combined_' datestr(now,'yyyy-mm-dd') '.fig'])



% Paper plots
% Power plots
figure('Position', [100 100 1100 850])
tiles = tiledlayout(3, 3,'TileSpacing','compact','Padding','compact');

% All Error 
nexttile(1)
contourf(data.time,data.freq,squeeze(Power_E_All_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Error'; 'Frequency (Hz)'}, 'FontWeight', 'bold')
title('Arrow')

nexttile(2)
contourf(data.time,data.freq,squeeze(Power_E_All_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
title('Unpleasant')

nexttile(3)
contourf(data.time,data.freq,squeeze(Power_E_All_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
title('Pleasant')
colorbar('east')

% All Correct 
nexttile(4)
contourf(data.time,data.freq,squeeze(Power_C_All_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Correct'; 'Frequency (Hz)'}, 'FontWeight', 'bold')

nexttile(5)
contourf(data.time,data.freq,squeeze(Power_C_All_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])

nexttile(6)
contourf(data.time,data.freq,squeeze(Power_C_All_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
colorbar('east')

% All Differences 
nexttile(7)
contourf(data.time,data.freq,squeeze(Power_E_All_Arr(chan_n_tf_theta,:,:) - Power_C_All_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Difference'; 'Frequency (Hz)'}, 'FontWeight', 'bold')
xlabel('Time (ms)', 'FontWeight', 'bold')

nexttile(8)
contourf(data.time,data.freq,squeeze(Power_E_All_Neg(chan_n_tf_theta,:,:) - Power_C_All_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
xlabel('Time (ms)', 'FontWeight', 'bold')

nexttile(9)
contourf(data.time,data.freq,squeeze(Power_E_All_Pos(chan_n_tf_theta,:,:) - Power_C_All_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
xlabel('Time (ms)', 'FontWeight', 'bold')
colorbar('east')

title(tiles, 'Power (dB)', 'FontSize', 16, 'FontWeight', 'Bold')

saveas(gcf, [fPath filesep 'Resp_Theta_Power_All_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_Theta_Power_All_' datestr(now,'yyyy-mm-dd') '.fig'])


% ITPC plots
figure('Position', [100 100 1100 850])
tiles = tiledlayout(3, 3,'TileSpacing','compact','Padding','compact');

% All Error 
nexttile(1)
contourf(data.time,data.freq,squeeze(ITPC_E_All_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Error'; 'Frequency (Hz)'}, 'FontWeight', 'bold')
title('Arrow')

nexttile(2)
contourf(data.time,data.freq,squeeze(ITPC_E_All_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
title('Unpleasant')

nexttile(3)
contourf(data.time,data.freq,squeeze(ITPC_E_All_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
title('Pleasant')
colorbar('east')

% All Correct 
nexttile(4)
contourf(data.time,data.freq,squeeze(ITPC_C_All_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Correct'; 'Frequency (Hz)'}, 'FontWeight', 'bold')

nexttile(5)
contourf(data.time,data.freq,squeeze(ITPC_C_All_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])

nexttile(6)
contourf(data.time,data.freq,squeeze(ITPC_C_All_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
colorbar('east')

% All Differences 
nexttile(7)
contourf(data.time,data.freq,squeeze(ITPC_E_All_Arr(chan_n_tf_theta,:,:) - ITPC_C_All_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-.25 .25],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Difference'; 'Frequency (Hz)'}, 'FontWeight', 'bold')
xlabel('Time (ms)', 'FontWeight', 'bold')

nexttile(8)
contourf(data.time,data.freq,squeeze(ITPC_E_All_Neg(chan_n_tf_theta,:,:) - ITPC_C_All_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-.25 .25],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
xlabel('Time (ms)', 'FontWeight', 'bold')

nexttile(9)
contourf(data.time,data.freq,squeeze(ITPC_E_All_Pos(chan_n_tf_theta,:,:) - ITPC_C_All_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-.25 .25],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
xlabel('Time (ms)', 'FontWeight', 'bold')
colorbar('east')

title(tiles, 'ITPC', 'FontSize', 16, 'FontWeight', 'Bold')

saveas(gcf, [fPath filesep 'Resp_Theta_ITPC_All_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_Theta_ITPC_All_' datestr(now,'yyyy-mm-dd') '.fig'])


% Power topography
figure('Position', [100 100 1100 850])

limit_Arr = max(abs(mean(Power_E_All_Arr(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')));
limit_Neg = max(abs(mean(Power_E_All_Neg(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')));
limit_Pos = max(abs(mean(Power_E_All_Pos(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')));
limit = max([limit_Arr, limit_Neg, limit_Pos]);
topo_limits = [-limit limit];

% All Error 
subplot(3,3,1)
topoplot((mean(Power_E_All_Arr(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')), data.chan_tf, 'maplimits', topo_limits, 'plotrad',.53);
text(-.8, 0, 'Error', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 14)
colormap(jet)
title('Arrow', 'FontSize', 14)

subplot(3,3,2)
topoplot((mean(Power_E_All_Neg(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')), data.chan_tf, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)
title('Unpleasant', 'FontSize', 14)

subplot(3,3,3)
topoplot((mean(Power_E_All_Pos(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')), data.chan_tf, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)
title('Pleasant', 'FontSize', 14)

c = colorbar('eastoutside');
pos = get(gca, 'Position');
c.Position = [pos(1) + pos(3) + .01, pos(2), .02, pos(4)];

% All Correct 
subplot(3,3,4)
topoplot((mean(Power_C_All_Arr(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')), data.chan_tf, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)
text(-.8, 0, 'Correct', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 14)

subplot(3,3,5)
topoplot((mean(Power_C_All_Neg(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')), data.chan_tf, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

subplot(3,3,6)
topoplot((mean(Power_C_All_Pos(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')), data.chan_tf, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

c = colorbar('eastoutside');
pos = get(gca, 'Position');
c.Position = [pos(1) + pos(3) + .01, pos(2), .02, pos(4)];

% All Differences 
subplot(3,3,7)
topoplot((mean(Power_E_All_Arr(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')- mean(Power_C_All_Arr (:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')),data.chan_tf, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)
text(-.8, 0, 'Difference', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 14)

subplot(3,3,8)
topoplot((mean(Power_E_All_Neg(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan') - mean(Power_C_All_Neg (:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')),data.chan_tf, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

subplot(3,3,9)
topoplot((mean(Power_E_All_Pos(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan') - mean(Power_C_All_Pos (:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')),data.chan_tf, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

c = colorbar('eastoutside');
pos = get(gca, 'Position');
c.Position = [pos(1) + pos(3) + .01, pos(2), .02, pos(4)];

annotation('textbox', [.52, 1, 0, 0], 'String', 'Theta Power (dB)', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'Bold', 'FitBoxToText', 'on', 'EdgeColor', 'none')

saveas(gcf, [fPath filesep 'Resp_Theta_Power_Topography_All_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_Theta_Power_Topography_All_' datestr(now,'yyyy-mm-dd') '.fig'])


% ITPC topography
figure('Position', [100 100 1100 850])

limit_Arr = max(abs(mean(ITPC_E_All_Arr(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')));
limit_Neg = max(abs(mean(ITPC_E_All_Neg(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')));
limit_Pos = max(abs(mean(ITPC_E_All_Pos(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')));
limit = max([limit_Arr, limit_Neg, limit_Pos]);
topo_limits = [0 limit];

% All Error 
subplot(3,3,1)
topoplot((mean(ITPC_E_All_Arr(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')), data.chan_tf, 'maplimits', topo_limits, 'plotrad',.53);
text(-.8, 0, 'Error', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 14)
colormap(jet)
title('Arrow', 'FontSize', 14)

subplot(3,3,2)
topoplot((mean(ITPC_E_All_Neg(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')), data.chan_tf, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)
title('Unpleasant', 'FontSize', 14)

subplot(3,3,3)
topoplot((mean(ITPC_E_All_Pos(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')), data.chan_tf, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)
title('Pleasant', 'FontSize', 14)

c = colorbar('eastoutside');
pos = get(gca, 'Position');
c.Position = [pos(1) + pos(3) + .01, pos(2), .02, pos(4)];

% All Correct 
subplot(3,3,4)
topoplot((mean(ITPC_C_All_Arr(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')), data.chan_tf, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)
text(-.8, 0, 'Correct', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 14)

subplot(3,3,5)
topoplot((mean(ITPC_C_All_Neg(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')), data.chan_tf, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

subplot(3,3,6)
topoplot((mean(ITPC_C_All_Pos(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')), data.chan_tf, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

c = colorbar('eastoutside');
pos = get(gca, 'Position');
c.Position = [pos(1) + pos(3) + .01, pos(2), .02, pos(4)];

% All Differences 
topo_limits = [-limit/4 limit/4];
subplot(3,3,7)
topoplot((mean(ITPC_E_All_Arr(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')- mean(ITPC_C_All_Arr (:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')),data.chan_tf, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)
text(-.8, 0, 'Difference', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 14)

subplot(3,3,8)
topoplot((mean(ITPC_E_All_Neg(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan') - mean(ITPC_C_All_Neg (:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')),data.chan_tf, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

subplot(3,3,9)
topoplot((mean(ITPC_E_All_Pos(:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan') - mean(ITPC_C_All_Pos (:, freq_tf_max_lower_idx: freq_tf_max_upper_idx, time_tf_max_lower_idx: time_tf_max_upper_idx), [2, 3], 'omitnan')),data.chan_tf, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

c = colorbar('eastoutside');
pos = get(gca, 'Position');
c.Position = [pos(1) + pos(3) + .01, pos(2), .02, pos(4)];

annotation('textbox', [.52, 1, 0, 0], 'String', 'Theta ITPC', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'Bold', 'FitBoxToText', 'on', 'EdgeColor', 'none')

saveas(gcf, [fPath filesep 'Resp_Theta_ITPC_Topography_All_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_Theta_ITPC_Topography_All_' datestr(now,'yyyy-mm-dd') '.fig'])


% ERN - matched
set(0,'defaultAxesFontSize',14, 'defaultTextFontSize',14)

figure('Position', [100 100 1100 600])
tiles = tiledlayout(2, 3,'TileSpacing','compact','Padding','compact');

% All Error and Correct 
nexttile(1)
plot(data.time,squeeze(ERP_E_All_Arr(chan_n_ern,:)), data.time, squeeze(ERP_C_All_Arr(chan_n_ern,:)))
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 12])
patch([ERN_lower ERN_lower ERN_upper ERN_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
ylabel({'Error and Correct'; 'Amplitude (\muV)'}, 'FontWeight', 'bold')
title('Arrow')

nexttile(2)
plot(data.time,squeeze(ERP_E_All_Neg(chan_n_ern,:)), data.time, squeeze(ERP_C_All_Neg(chan_n_ern,:)))
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 12])
patch([ERN_lower ERN_lower ERN_upper ERN_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
title('Unpleasant')

nexttile(3)
plot(data.time,squeeze(ERP_E_All_Pos(chan_n_ern,:)), data.time, squeeze(ERP_C_All_Pos(chan_n_ern,:)))
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 12])
patch([ERN_lower ERN_lower ERN_upper ERN_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
legend('Error', 'Correct', 'Location', 'southeast')
title('Pleasant')

% All Differences 
nexttile(4)
plot(data.time,squeeze(ERP_E_All_Arr(chan_n_ern,:)-ERP_C_All_Arr(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 12])
patch([ERN_lower ERN_lower ERN_upper ERN_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
ylabel({'Difference'; 'Amplitude (\muV)'}, 'FontWeight', 'bold')
xlabel('Time (ms)', 'FontWeight', 'bold')

nexttile(5)
plot(data.time,squeeze(ERP_E_All_Neg(chan_n_ern,:)-ERP_C_All_Neg(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 12])
patch([ERN_lower ERN_lower ERN_upper ERN_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
xlabel('Time (ms)', 'FontWeight', 'bold')

nexttile(6)
plot(data.time,squeeze(ERP_E_All_Pos(chan_n_ern,:)-ERP_C_All_Pos(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 12])
patch([ERN_lower ERN_lower ERN_upper ERN_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
legend('Difference')
xlabel('Time (ms)', 'FontWeight', 'bold')

title(tiles, 'ERN, CRN, and \DeltaERN', 'FontSize', 16, 'FontWeight', 'Bold')

saveas(gcf, [fPath filesep 'Resp_ERN_All_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_ERN_All_' datestr(now,'yyyy-mm-dd') '.fig'])



% ERN – Non-Matched
figure('Position', [100 100 1100 600])
tiles = tiledlayout(2, 3,'TileSpacing','compact','Padding','compact');

% All Error and Correct 
nexttile(1)
plot(data.time,squeeze(ERP_E_All_Arr(chan_n_ern,:)), data.time, squeeze(ERP_C_non_All_Arr(chan_n_ern,:)))
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 12])
patch([ERN_non_lower ERN_non_lower ERN_non_upper ERN_non_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
ylabel({'Error and Correct'; 'Amplitude (\muV)'}, 'FontWeight', 'bold')
title('Arrow')

nexttile(2)
plot(data.time,squeeze(ERP_E_All_Neg(chan_n_ern,:)), data.time, squeeze(ERP_C_non_All_Neg(chan_n_ern,:)))
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 12])
patch([ERN_non_lower ERN_non_lower ERN_non_upper ERN_non_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
title('Unpleasant')

nexttile(3)
plot(data.time,squeeze(ERP_E_All_Pos(chan_n_ern,:)), data.time, squeeze(ERP_C_non_All_Pos(chan_n_ern,:)))
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 12])
patch([ERN_non_lower ERN_non_lower ERN_non_upper ERN_non_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
legend('Error', 'Correct (All Trials)', 'Location', 'southeast')
title('Pleasant')

% All Differences 
nexttile(4)
plot(data.time,squeeze(ERP_E_All_Arr(chan_n_ern,:)-ERP_C_non_All_Arr(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 12])
patch([ERN_non_lower ERN_non_lower ERN_non_upper ERN_non_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
ylabel({'Difference'; 'Amplitude (\muV)'}, 'FontWeight', 'bold')

nexttile(5)
plot(data.time,squeeze(ERP_E_All_Neg(chan_n_ern,:)-ERP_C_non_All_Neg(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 12])
patch([ERN_non_lower ERN_non_lower ERN_non_upper ERN_non_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])

nexttile(6)
plot(data.time,squeeze(ERP_E_All_Pos(chan_n_ern,:)-ERP_C_non_All_Pos(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 12])
patch([ERN_non_lower ERN_non_lower ERN_non_upper ERN_non_upper], [-12 15 15 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
legend('Difference')

title(tiles, 'ERN, CRN (All Trials), and \DeltaERN', 'FontSize', 16, 'FontWeight', 'Bold')

saveas(gcf, [fPath filesep 'Resp_ERN_non_All_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_ERN_non_All_' datestr(now,'yyyy-mm-dd') '.fig'])


% ERN topography
figure('Position', [100 100 1100 850])

limit_Arr = max(abs(mean(ERP_C_All_Arr(:, ERN_lower_idx: ERN_upper_idx), 2, 'omitnan')));
limit_Neg = max(abs(mean(ERP_C_All_Neg(:, ERN_lower_idx: ERN_upper_idx), 2, 'omitnan')));
limit_Pos = max(abs(mean(ERP_C_All_Pos(:, ERN_lower_idx: ERN_upper_idx), 2, 'omitnan')));
limit = max([limit_Arr, limit_Neg, limit_Pos]);
topo_limits = [-limit limit];

% All Error 
subplot(3,3,1)
topoplot((mean(ERP_E_All_Arr(:, ERN_lower_idx: ERN_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
text(-.8, 0, 'Error', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 14)
colormap(jet)
title('Arrow', 'FontSize', 14)

subplot(3,3,2)
topoplot((mean(ERP_E_All_Neg(:, ERN_lower_idx: ERN_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)
title('Unpleasant', 'FontSize', 14)

subplot(3,3,3)
topoplot((mean(ERP_E_All_Pos(:, ERN_lower_idx: ERN_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)
title('Pleasant', 'FontSize', 14)

c = colorbar('eastoutside');
pos = get(gca, 'Position');
c.Position = [pos(1) + pos(3) + .01, pos(2), .02, pos(4)];

% All Correct 
subplot(3,3,4)
topoplot((mean(ERP_C_All_Arr(:, ERN_lower_idx: ERN_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
text(-.8, 0, 'Correct', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 14)
colormap(jet)

subplot(3,3,5)
topoplot((mean(ERP_C_All_Neg(:, ERN_lower_idx: ERN_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

subplot(3,3,6)
topoplot((mean(ERP_C_All_Pos(:, ERN_lower_idx: ERN_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

c = colorbar('eastoutside');
pos = get(gca, 'Position');
c.Position = [pos(1) + pos(3) + .01, pos(2), .02, pos(4)];

% All differences 
subplot(3,3,7)
topoplot((mean(ERP_E_All_Arr(:, ERN_lower_idx: ERN_upper_idx), 2, 'omitnan') - mean(ERP_C_All_Arr(:, ERN_lower_idx: ERN_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
text(-.8, 0, 'Difference', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 14)
colormap(jet)

subplot(3,3,8)
topoplot((mean(ERP_E_All_Neg(:, ERN_lower_idx: ERN_upper_idx), 2, 'omitnan') - mean(ERP_C_All_Neg(:, ERN_lower_idx: ERN_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

subplot(3,3,9)
topoplot((mean(ERP_E_All_Pos(:, ERN_lower_idx: ERN_upper_idx), 2, 'omitnan') - mean(ERP_C_All_Pos(:, ERN_lower_idx: ERN_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

c = colorbar('eastoutside');
pos = get(gca, 'Position');
c.Position = [pos(1) + pos(3) + .01, pos(2), .02, pos(4)];

annotation('textbox', [.52, 1, 0, 0], 'String', 'ERN, CRN, and \DeltaERN', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'Bold', 'FitBoxToText', 'on', 'EdgeColor', 'none')

saveas(gcf, [fPath filesep 'Resp_ERN_Topography_All_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_ERN_Topography_All_' datestr(now,'yyyy-mm-dd') '.fig'])



% ERN Non-Matched topography
figure('Position', [100 100 1100 850])

limit_Arr = max(abs(mean(ERP_C_non_All_Arr(:, ERN_non_lower_idx: ERN_non_upper_idx), 2, 'omitnan')));
limit_Neg = max(abs(mean(ERP_C_non_All_Neg(:, ERN_non_lower_idx: ERN_non_upper_idx), 2, 'omitnan')));
limit_Pos = max(abs(mean(ERP_C_non_All_Pos(:, ERN_non_lower_idx: ERN_non_upper_idx), 2, 'omitnan')));
limit = max([limit_Arr, limit_Neg, limit_Pos]);
topo_limits = [-limit limit];

% All Error 
subplot(3,3,1)
topoplot((mean(ERP_E_All_Arr(:, ERN_non_lower_idx: ERN_non_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
text(-.8, 0, 'Error', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 14)
colormap(jet)
title('Arrow', 'FontSize', 14)

subplot(3,3,2)
topoplot((mean(ERP_E_All_Neg(:, ERN_non_lower_idx: ERN_non_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)
title('Unpleasant', 'FontSize', 14)

subplot(3,3,3)
topoplot((mean(ERP_E_All_Pos(:, ERN_non_lower_idx: ERN_non_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)
title('Pleasant', 'FontSize', 14)

c = colorbar('eastoutside');
pos = get(gca, 'Position');
c.Position = [pos(1) + pos(3) + .01, pos(2), .02, pos(4)];

% All Correct 
subplot(3,3,4)
topoplot((mean(ERP_C_non_All_Arr(:, ERN_non_lower_idx: ERN_non_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
text(-.8, 0, 'Correct', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 14)
colormap(jet)

subplot(3,3,5)
topoplot((mean(ERP_C_non_All_Neg(:, ERN_non_lower_idx: ERN_non_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

subplot(3,3,6)
topoplot((mean(ERP_C_non_All_Pos(:, ERN_non_lower_idx: ERN_non_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

c = colorbar('eastoutside');
pos = get(gca, 'Position');
c.Position = [pos(1) + pos(3) + .01, pos(2), .02, pos(4)];

% All Differences 
subplot(3,3,7)
topoplot((mean(ERP_E_All_Arr(:, ERN_non_lower_idx: ERN_non_upper_idx), 2, 'omitnan') - mean(ERP_C_non_All_Arr(:, ERN_non_lower_idx: ERN_non_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
text(-.8, 0, 'Difference', 'HorizontalAlignment', 'center', 'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 14)
colormap(jet)

subplot(3,3,8)
topoplot((mean(ERP_E_All_Neg(:, ERN_non_lower_idx: ERN_non_upper_idx), 2, 'omitnan') - mean(ERP_C_non_All_Neg(:, ERN_non_lower_idx: ERN_non_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

subplot(3,3,9)
topoplot((mean(ERP_E_All_Pos(:, ERN_non_lower_idx: ERN_non_upper_idx), 2, 'omitnan') - mean(ERP_C_non_All_Pos(:, ERN_non_lower_idx: ERN_non_upper_idx), 2, 'omitnan')), data.chan_erp, 'maplimits', topo_limits, 'plotrad',.53);
colormap(jet)

c = colorbar('eastoutside');
pos = get(gca, 'Position');
c.Position = [pos(1) + pos(3) + .01, pos(2), .02, pos(4)];

annotation('textbox', [.52, 1, 0, 0], 'String', 'ERN, CRN (All Trials), and \DeltaERN', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'Bold', 'FitBoxToText', 'on', 'EdgeColor', 'none')

saveas(gcf, [fPath filesep 'Resp_Theta_ERP_non_Topography_All_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_Theta_ERP_non_Topography_All_' datestr(now,'yyyy-mm-dd') '.fig'])


%% Group Differences
% Theta
% Arrow Power
set(0,'defaultAxesFontSize',14, 'defaultTextFontSize',14)

figure('Position', [100 100 1100 850])
tiles = tiledlayout(3, 3,'TileSpacing','compact','Padding','compact');

% Arrow Error by Group
nexttile(1)
contourf(data.time,data.freq,squeeze(Power_E_HC_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Error'; 'Frequency (Hz)'}, 'FontWeight', 'bold')
title('HC')

nexttile(2)
contourf(data.time,data.freq,squeeze(Power_E_BD_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
title('BD')

nexttile(3)
contourf(data.time,data.freq,squeeze(Power_E_SZ_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
title('SZ')


colorbar('east')

% Arrow Correct by Group
nexttile(4)
contourf(data.time,data.freq,squeeze(Power_C_HC_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Correct'; 'Frequency (Hz)'}, 'FontWeight', 'bold')

nexttile(5)
contourf(data.time,data.freq,squeeze(Power_C_BD_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])

nexttile(6)
contourf(data.time,data.freq,squeeze(Power_C_SZ_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])

colorbar('east')

% Arrow differences by Group
nexttile(7)
contourf(data.time,data.freq,squeeze(Power_E_HC_Arr(chan_n_tf_theta,:,:) - Power_C_HC_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Difference'; 'Frequency (Hz)'}, 'FontWeight', 'bold')
xlabel('Time (ms)', 'FontWeight', 'bold')

nexttile(8)
contourf(data.time,data.freq,squeeze(Power_E_BD_Arr(chan_n_tf_theta,:,:) - Power_C_BD_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
xlabel('Time (ms)', 'FontWeight', 'bold')

nexttile(9)
contourf(data.time,data.freq,squeeze(Power_E_SZ_Arr(chan_n_tf_theta,:,:) - Power_C_SZ_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
xlabel('Time (ms)', 'FontWeight', 'bold')

colorbar('east')

title(tiles, 'Arrow, Power (dB)', 'FontSize', 16, 'FontWeight', 'Bold')

saveas(gcf, [fPath filesep 'Resp_Theta_Power_Groups_Arr_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_Theta_Power_Groups_Arr_' datestr(now,'yyyy-mm-dd') '.fig'])


% Unpleasant Power
figure('Position', [100 100 1100 850])
tiles = tiledlayout(3, 3,'TileSpacing','compact','Padding','compact');

% Unpleasant Error by Group
nexttile(1)
contourf(data.time,data.freq,squeeze(Power_E_HC_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Error'; 'Frequency (Hz)'}, 'FontWeight', 'bold')
title('HC')

nexttile(2)
contourf(data.time,data.freq,squeeze(Power_E_BD_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
title('BD')

nexttile(3)
contourf(data.time,data.freq,squeeze(Power_E_SZ_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
title('SZ')

colorbar('east')

% Unpleasant Correct by Group
nexttile(4)
contourf(data.time,data.freq,squeeze(Power_C_HC_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Correct'; 'Frequency (Hz)'}, 'FontWeight', 'bold')

nexttile(5)
contourf(data.time,data.freq,squeeze(Power_C_BD_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])

nexttile(6)
contourf(data.time,data.freq,squeeze(Power_C_SZ_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])

colorbar('east')

% Unpleasant Differences by Group
nexttile(7)
contourf(data.time,data.freq,squeeze(Power_E_HC_Neg(chan_n_tf_theta,:,:) - Power_C_HC_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Difference'; 'Frequency (Hz)'}, 'FontWeight', 'bold')
xlabel('Time (ms)', 'FontWeight', 'bold')

nexttile(8)
contourf(data.time,data.freq,squeeze(Power_E_BD_Neg(chan_n_tf_theta,:,:) - Power_C_BD_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
xlabel('Time (ms)', 'FontWeight', 'bold')

nexttile(9)
contourf(data.time,data.freq,squeeze(Power_E_SZ_Neg(chan_n_tf_theta,:,:) - Power_C_SZ_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
xlabel('Time (ms)', 'FontWeight', 'bold')

colorbar('east')

title(tiles, 'Unpleasant, Power (dB)', 'FontSize', 16, 'FontWeight', 'Bold')

saveas(gcf, [fPath filesep 'Resp_Theta_Power_Groups_Neg_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_Theta_Power_Groups_Neg_' datestr(now,'yyyy-mm-dd') '.fig'])


% Pleasant Power
figure('Position', [100 100 1100 850])
tiles = tiledlayout(3, 3,'TileSpacing','compact','Padding','compact');

% Pleasant Error by Group
nexttile(1)
contourf(data.time,data.freq,squeeze(Power_E_HC_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Error'; 'Frequency (Hz)'}, 'FontWeight', 'bold')
title('HC')

nexttile(2)
contourf(data.time,data.freq,squeeze(Power_E_BD_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
title('BD')

nexttile(3)
contourf(data.time,data.freq,squeeze(Power_E_SZ_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
title('SZ')

colorbar('east')


% Correct Error by Group
nexttile(4)
contourf(data.time,data.freq,squeeze(Power_C_HC_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Correct'; 'Frequency (Hz)'}, 'FontWeight', 'bold')

nexttile(5)
contourf(data.time,data.freq,squeeze(Power_C_BD_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])

nexttile(6)
contourf(data.time,data.freq,squeeze(Power_C_SZ_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])

colorbar('east')

% Pleasant Difference by Group
nexttile(7)
contourf(data.time,data.freq,squeeze(Power_E_HC_Pos(chan_n_tf_theta,:,:) - Power_C_HC_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Difference'; 'Frequency (Hz)'}, 'FontWeight', 'bold')
xlabel('Time (ms)', 'FontWeight', 'bold')

nexttile(8)
contourf(data.time,data.freq,squeeze(Power_E_BD_Pos(chan_n_tf_theta,:,:) - Power_C_BD_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
xlabel('Time (ms)', 'FontWeight', 'bold')

nexttile(9)
contourf(data.time,data.freq,squeeze(Power_E_SZ_Pos(chan_n_tf_theta,:,:) - Power_C_SZ_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-4 4],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
xlabel('Time (ms)', 'FontWeight', 'bold')

colorbar('east')

title(tiles, 'Pleasant, Power (dB)', 'FontSize', 16, 'FontWeight', 'Bold')

saveas(gcf, [fPath filesep 'Resp_Theta_Power_Groups_Pos_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_Theta_Power_Groups_Pos_' datestr(now,'yyyy-mm-dd') '.fig'])




% Arrow ITPC
figure('Position', [100 100 1100 850])
tiles = tiledlayout(3, 3,'TileSpacing','compact','Padding','compact');

% Arrow Error by Group
nexttile(1)
contourf(data.time,data.freq,squeeze(ITPC_E_HC_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Error'; 'Frequency (Hz)'}, 'FontWeight', 'bold')
title('HC')

nexttile(2)
contourf(data.time,data.freq,squeeze(ITPC_E_BD_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
title('BD')

nexttile(3)
contourf(data.time,data.freq,squeeze(ITPC_E_SZ_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
title('SZ')

colorbar('east')


% Arrow Correct by Group
nexttile(4)
contourf(data.time,data.freq,squeeze(ITPC_C_HC_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Correct'; 'Frequency (Hz)'}, 'FontWeight', 'bold')

nexttile(5)
contourf(data.time,data.freq,squeeze(ITPC_C_BD_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])

nexttile(6)
contourf(data.time,data.freq,squeeze(ITPC_C_SZ_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])

colorbar('east')

% Arrow Differences by Group
nexttile(7)
contourf(data.time,data.freq,squeeze(ITPC_E_HC_Arr(chan_n_tf_theta,:,:) - ITPC_C_HC_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-.25 .25],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Difference'; 'Frequency (Hz)'}, 'FontWeight', 'bold')
xlabel('Time (ms)', 'FontWeight', 'bold')

nexttile(8)
contourf(data.time,data.freq,squeeze(ITPC_E_BD_Arr(chan_n_tf_theta,:,:) - ITPC_C_BD_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-.25 .25],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
xlabel('Time (ms)', 'FontWeight', 'bold')

nexttile(9)
contourf(data.time,data.freq,squeeze(ITPC_E_SZ_Arr(chan_n_tf_theta,:,:) - ITPC_C_SZ_Arr(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-.25 .25],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
xlabel('Time (ms)', 'FontWeight', 'bold')

colorbar('east')

title(tiles, 'Arrow, ITPC', 'FontSize', 16, 'FontWeight', 'Bold')

saveas(gcf, [fPath filesep 'Resp_Theta_ITPC_Groups_Arr_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_Theta_ITPC_Groups_Arr_' datestr(now,'yyyy-mm-dd') '.fig'])

% Arrow ITPC
figure('Position', [100 100 1100 850])
tiles = tiledlayout(3, 3,'TileSpacing','compact','Padding','compact');

% Unpleasant Error by Group
nexttile(1)
contourf(data.time,data.freq,squeeze(ITPC_E_HC_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Error'; 'Frequency (Hz)'}, 'FontWeight', 'bold')
title('HC')

nexttile(2)
contourf(data.time,data.freq,squeeze(ITPC_E_BD_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
title('BD')

nexttile(3)
contourf(data.time,data.freq,squeeze(ITPC_E_SZ_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
title('SZ')

colorbar('east')


% Unpleasant Correct by Group
nexttile(4)
contourf(data.time,data.freq,squeeze(ITPC_C_HC_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Correct'; 'Frequency (Hz)'}, 'FontWeight', 'bold')

nexttile(5)
contourf(data.time,data.freq,squeeze(ITPC_C_BD_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])

nexttile(6)
contourf(data.time,data.freq,squeeze(ITPC_C_SZ_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])

colorbar('east')

% Unpleasant Differences by Group
nexttile(7)
contourf(data.time,data.freq,squeeze(ITPC_E_HC_Neg(chan_n_tf_theta,:,:) - ITPC_C_HC_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-.25 .25],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Difference'; 'Frequency (Hz)'}, 'FontWeight', 'bold')
xlabel('Time (ms)', 'FontWeight', 'bold')

nexttile(8)
contourf(data.time,data.freq,squeeze(ITPC_E_BD_Neg(chan_n_tf_theta,:,:) - ITPC_C_BD_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-.25 .25],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
xlabel('Time (ms)', 'FontWeight', 'bold')

nexttile(9)
contourf(data.time,data.freq,squeeze(ITPC_E_SZ_Neg(chan_n_tf_theta,:,:) - ITPC_C_SZ_Neg(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-.25 .25],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
xlabel('Time (ms)', 'FontWeight', 'bold')

colorbar('east')

title(tiles, 'Unpleasant, ITPC', 'FontSize', 16, 'FontWeight', 'Bold')

saveas(gcf, [fPath filesep 'Resp_Theta_ITPC_Groups_Neg_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_Theta_ITPC_Groups_Neg_' datestr(now,'yyyy-mm-dd') '.fig'])



% Pleasant by Group
figure('Position', [100 100 1100 850])
tiles = tiledlayout(3, 3,'TileSpacing','compact','Padding','compact');

% Pleasant Error by Group
nexttile(1)
contourf(data.time,data.freq,squeeze(ITPC_E_HC_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Error'; 'Frequency (Hz)'}, 'FontWeight', 'bold')
title('HC')

nexttile(2)
contourf(data.time,data.freq,squeeze(ITPC_E_BD_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
title('BD')

nexttile(3)
contourf(data.time,data.freq,squeeze(ITPC_E_SZ_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
title('SZ')

colorbar('east')


% Pleasant Correct by Group
nexttile(4)
contourf(data.time,data.freq,squeeze(ITPC_C_HC_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Correct'; 'Frequency (Hz)'}, 'FontWeight', 'bold')

nexttile(5)
contourf(data.time,data.freq,squeeze(ITPC_C_BD_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])

nexttile(6)
contourf(data.time,data.freq,squeeze(ITPC_C_SZ_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[0 .5],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])

colorbar('east')

% Positivwe Differences by Group
nexttile(7)
contourf(data.time,data.freq,squeeze(ITPC_E_HC_Pos(chan_n_tf_theta,:,:) - ITPC_C_HC_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-.25 .25],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
ylabel({'Difference'; 'Frequency (Hz)'}, 'FontWeight', 'bold')
xlabel('Time (ms)', 'FontWeight', 'bold')

nexttile(8)
contourf(data.time,data.freq,squeeze(ITPC_E_BD_Pos(chan_n_tf_theta,:,:) - ITPC_C_BD_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-.25 .25],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
xlabel('Time (ms)', 'FontWeight', 'bold')

nexttile(9)
contourf(data.time,data.freq,squeeze(ITPC_E_SZ_Pos(chan_n_tf_theta,:,:) - ITPC_C_SZ_Pos(chan_n_tf_theta,:,:)),100,'linecolor','none')
set(gca,'clim',[-.25 .25],'xlim',[-500 1000], 'ylim', [min_freq max_freq], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colormap(jet)
rectangle('position', [time_tf_max_lower freq_tf_max_lower time_tf_max_width freq_tf_max_width])
xlabel('Time (ms)', 'FontWeight', 'bold')

colorbar('east')

title(tiles, 'Pleasant, ITPC', 'FontSize', 16, 'FontWeight', 'Bold')

saveas(gcf, [fPath filesep 'Resp_Theta_ITPC_Groups_Pos_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_Theta_ITPC_Groups_Pos_' datestr(now,'yyyy-mm-dd') '.fig'])


% ERN - matched by Group
figure('Position', [100 100 1100 600])
tiles = tiledlayout(3, 3,'TileSpacing','compact','Padding','compact');

% Arrow by Group
nexttile(1)
plot(data.time,squeeze(ERP_E_HC_Arr(chan_n_ern,:)), '--', data.time, squeeze(ERP_C_HC_Arr(chan_n_ern,:)), '-.', data.time, squeeze(ERP_E_HC_Arr(chan_n_ern,:)-ERP_C_HC_Arr(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 18])
patch([ERN_lower ERN_lower ERN_upper ERN_upper], [-12 18 18 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
ylabel({'Arrow'; 'Amplitude (\muV)'}, 'FontWeight', 'bold')
title('HC')

nexttile(2)
plot(data.time,squeeze(ERP_E_BD_Arr(chan_n_ern,:)), '--', data.time, squeeze(ERP_C_BD_Arr(chan_n_ern,:)), '-.', data.time, squeeze(ERP_E_BD_Arr(chan_n_ern,:)-ERP_C_BD_Arr(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim',[-12, 18])
patch([ERN_lower ERN_lower ERN_upper ERN_upper], [-12 18 18 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
title('BD')

nexttile(3)
plot(data.time,squeeze(ERP_E_SZ_Arr(chan_n_ern,:)), '--', data.time, squeeze(ERP_C_SZ_Arr(chan_n_ern,:)), '-.', data.time, squeeze(ERP_E_SZ_Arr(chan_n_ern,:)-ERP_C_SZ_Arr(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 18])
patch([ERN_lower ERN_lower ERN_upper ERN_upper], [-12 18 18 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
title('SZ')
legend('Error', 'Correct', 'Difference', 'Location', 'southeast', 'Fontsize', 8)

% Unpleasant by Group
nexttile(4)
plot(data.time,squeeze(ERP_E_HC_Neg(chan_n_ern,:)), '--', data.time, squeeze(ERP_C_HC_Neg(chan_n_ern,:)), '-.', data.time, squeeze(ERP_E_HC_Neg(chan_n_ern,:)-ERP_C_HC_Neg(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 18])
patch([ERN_lower ERN_lower ERN_upper ERN_upper], [-12 18 18 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
ylabel({'Unpleasant'; 'Amplitude (\muV)'}, 'FontWeight', 'bold')

nexttile(5)
plot(data.time,squeeze(ERP_E_BD_Neg(chan_n_ern,:)), '--', data.time, squeeze(ERP_C_BD_Neg(chan_n_ern,:)), '-.', data.time, squeeze(ERP_E_BD_Neg(chan_n_ern,:)-ERP_C_BD_Neg(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 18])
patch([ERN_lower ERN_lower ERN_upper ERN_upper], [-12 18 18 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])

nexttile(6)
plot(data.time,squeeze(ERP_E_SZ_Neg(chan_n_ern,:)), '--', data.time, squeeze(ERP_C_SZ_Neg(chan_n_ern,:)), '-.', data.time, squeeze(ERP_E_SZ_Neg(chan_n_ern,:)-ERP_C_SZ_Neg(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 18])
patch([ERN_lower ERN_lower ERN_upper ERN_upper], [-12 18 18 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
legend('Error', 'Correct', 'Difference', 'Location', 'southeast', 'Fontsize', 8)

% Pleasant by Group
nexttile(7)
plot(data.time,squeeze(ERP_E_HC_Pos(chan_n_ern,:)), '--', data.time, squeeze(ERP_C_HC_Pos(chan_n_ern,:)), '-.', data.time, squeeze(ERP_E_HC_Pos(chan_n_ern,:)-ERP_C_HC_Pos(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 18])
patch([ERN_lower ERN_lower ERN_upper ERN_upper], [-12 18 18 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
ylabel({'Pleasant'; 'Amplitude (\muV)'}, 'FontWeight', 'bold')

nexttile(8)
plot(data.time,squeeze(ERP_E_BD_Pos(chan_n_ern,:)), '--', data.time, squeeze(ERP_C_BD_Pos(chan_n_ern,:)), '-.', data.time, squeeze(ERP_E_BD_Pos(chan_n_ern,:)-ERP_C_BD_Pos(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 18])
patch([ERN_lower ERN_lower ERN_upper ERN_upper], [-12 18 18 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])

nexttile(9)
plot(data.time,squeeze(ERP_E_SZ_Pos(chan_n_ern,:)), '--', data.time, squeeze(ERP_C_SZ_Pos(chan_n_ern,:)), '-.', data.time, squeeze(ERP_E_SZ_Pos(chan_n_ern,:)-ERP_C_SZ_Pos(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 18])
patch([ERN_lower ERN_lower ERN_upper ERN_upper], [-12 18 18 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
legend('Error', 'Correct', 'Difference', 'Location', 'southeast', 'Fontsize', 8)

title(tiles, 'ERN, CRN, and \DeltaERN', 'FontSize', 16, 'FontWeight', 'Bold')

saveas(gcf, [fPath filesep 'Resp_ERN_Groups_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_ERN_Groups_' datestr(now,'yyyy-mm-dd') '.fig'])

% ERN – non-matched by Group
figure('Position', [100 100 1100 600])
tiles = tiledlayout(3, 3,'TileSpacing','compact','Padding','compact');

% Arrow by Group
nexttile(1)
plot(data.time,squeeze(ERP_E_HC_Arr(chan_n_ern,:)), '--', data.time, squeeze(ERP_C_non_HC_Arr(chan_n_ern,:)), '-.', data.time, squeeze(ERP_E_HC_Arr(chan_n_ern,:)-ERP_C_non_HC_Arr(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 18])
patch([ERN_non_lower ERN_non_lower ERN_non_upper ERN_non_upper], [-12 18 18 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
ylabel({'Arrow'; 'Amplitude (\muV)'}, 'FontWeight', 'bold')
title('HC')

nexttile(2)
plot(data.time,squeeze(ERP_E_BD_Arr(chan_n_ern,:)), '--', data.time, squeeze(ERP_C_non_BD_Arr(chan_n_ern,:)), '-.', data.time, squeeze(ERP_E_BD_Arr(chan_n_ern,:)-ERP_C_non_BD_Arr(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 18])
patch([ERN_non_lower ERN_non_lower ERN_non_upper ERN_non_upper], [-12 18 18 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
title('BD')

nexttile(3)
plot(data.time,squeeze(ERP_E_SZ_Arr(chan_n_ern,:)), '--', data.time, squeeze(ERP_C_non_SZ_Arr(chan_n_ern,:)), '-.', data.time, squeeze(ERP_E_SZ_Arr(chan_n_ern,:)-ERP_C_non_SZ_Arr(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 18])
patch([ERN_non_lower ERN_non_lower ERN_non_upper ERN_non_upper], [-12 18 18 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
title('SZ')
legend('Error', 'Correct', 'Difference', 'Location', 'southeast', 'Fontsize', 8)

% Unpleasant by Group
nexttile(4)
plot(data.time,squeeze(ERP_E_HC_Neg(chan_n_ern,:)), '--', data.time, squeeze(ERP_C_non_HC_Neg(chan_n_ern,:)), '-.', data.time, squeeze(ERP_E_HC_Neg(chan_n_ern,:)-ERP_C_non_HC_Neg(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 18])
patch([ERN_non_lower ERN_non_lower ERN_non_upper ERN_non_upper], [-12 18 18 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
ylabel({'Unpleasant'; 'Amplitude (\muV)'}, 'FontWeight', 'bold')

nexttile(5)
plot(data.time,squeeze(ERP_E_BD_Neg(chan_n_ern,:)), '--', data.time, squeeze(ERP_C_non_BD_Neg(chan_n_ern,:)), '-.', data.time, squeeze(ERP_E_BD_Neg(chan_n_ern,:)-ERP_C_non_BD_Neg(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 18])
patch([ERN_non_lower ERN_non_lower ERN_non_upper ERN_non_upper], [-12 18 18 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])

nexttile(6)
plot(data.time,squeeze(ERP_E_SZ_Neg(chan_n_ern,:)), '--', data.time, squeeze(ERP_C_non_SZ_Neg(chan_n_ern,:)), '-.', data.time, squeeze(ERP_E_SZ_Neg(chan_n_ern,:)-ERP_C_non_SZ_Neg(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 18])
patch([ERN_non_lower ERN_non_lower ERN_non_upper ERN_non_upper], [-12 18 18 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
legend('Error', 'Correct', 'Difference', 'Location', 'southeast', 'Fontsize', 8)

% Pleasant by Group
nexttile(7)
plot(data.time,squeeze(ERP_E_HC_Pos(chan_n_ern,:)), '--', data.time, squeeze(ERP_C_non_HC_Pos(chan_n_ern,:)), '-.', data.time, squeeze(ERP_E_HC_Pos(chan_n_ern,:)-ERP_C_non_HC_Pos(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 18])
patch([ERN_non_lower ERN_non_lower ERN_non_upper ERN_non_upper], [-12 18 18 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
ylabel({'Pleasant'; 'Amplitude (\muV)'}, 'FontWeight', 'bold')

nexttile(8)
plot(data.time,squeeze(ERP_E_BD_Pos(chan_n_ern,:)), '--', data.time, squeeze(ERP_C_non_BD_Pos(chan_n_ern,:)), '-.', data.time, squeeze(ERP_E_BD_Pos(chan_n_ern,:)-ERP_C_non_BD_Pos(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 18])
patch([ERN_non_lower ERN_non_lower ERN_non_upper ERN_non_upper], [-12 18 18 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])

nexttile(9)
plot(data.time,squeeze(ERP_E_SZ_Pos(chan_n_ern,:)), '--', data.time, squeeze(ERP_C_non_SZ_Pos(chan_n_ern,:)), '-.', data.time, squeeze(ERP_E_SZ_Pos(chan_n_ern,:)-ERP_C_non_SZ_Pos(chan_n_ern,:)), "k")
set(gca, 'xlim', [-200, 600], 'ylim', [-12, 18])
patch([ERN_non_lower ERN_non_lower ERN_non_upper ERN_non_upper], [-12 18 18 -12], [.5 .5 .5], 'EdgeColor', 'none', 'Facealpha', .3)
yline(0, '-', 'Color', [.75 .75 .75])
legend('Error', 'Correct', 'Difference', 'Location', 'southeast', 'Fontsize', 8)

title(tiles, 'ERN, CRN (All Trials), and \DeltaERN', 'FontSize', 16, 'FontWeight', 'Bold')

saveas(gcf, [fPath filesep 'Resp_ERN_non_Groups_' datestr(now,'yyyy-mm-dd') '.tiff'])
savefig(gcf, [fPath filesep 'Resp_ERN_non_Groups_' datestr(now,'yyyy-mm-dd') '.fig'])


fprintf('Step 4 figure and data extraction are done');

