% Main File for Time-Frequency Analysis of Continuous EEG Signals
% Author: Ehsan Eqlimi, WAVES, UGent, Belgium
% Original Date: November 2020
% Last Updated: 15/12/2024 by EE

clc;
clear;
%close all;
% eeglab
%% Initialization Section
DoFaster=0;% Flag for faster execution
EEGDataPath='E:\Ehsan\Data_PP_Corrected\';% Path to the EEG data directory
addpath('E:\Ehsan\GitHub\ComplexPCA-PD');% Add custom analysis toolbox to the path

%% Add Required Toolboxes to Path
% Add the FieldTrip toolbox for EEG analysis and initialize its defaults
addpath('E:\Ehsan\RippleServerFiles\Toolbox\fieldtrip-20220729');
ft_defaults;

% Specify the domain for data selection
Domain='HC*avgall';% Options: 'HC*avgall' or 'HC*avgmast'
EDFDir=dir(fullfile(EEGDataPath,['*' Domain '.edf']));% List of EDF files in the specified domain
% Initialize an empty table for categorical data (if needed)
CatTable=[];

%% Epoching and Windowing Parameters
TimeRange=[0 6];% Time window for epoching in seconds
FreqRange=[7 13];
% BaselineWin=[0.3 0.6];% Baseline window in seconds
BaselineWin=[0 0.5];% Baseline window in seconds


% Time ranges for induced power analysis
TimeRange_PreSent_Ind=[1.5 2.0];% Previous selection: [1.5 2]
TimeRange_DuringSent_Ind=[2.5 3.5];% Previous selection: [2.5 4.5]
TimeRange_PostSent_Ind=[4.5 5.5];% Induced power analysis window
BaseLineNormtype='db';% Baseline normalization type (e.g., 'db')
% Define the time window range and step size
TimeWindowing=TimeRange(1):0.5:TimeRange(2); % Time range from 0 to 6 seconds with 0.5-second step

% Create pairs of start and end times
TimeWindowMatrix=[TimeWindowing(1:end-1)' TimeWindowing(2:end)'];
%% Time-Frequency Analysis Parameters
TFParam.pad=[];% Zero-padding (if any)
TFParam.keeptrials='yes';% Keep individual trials during analysis
TFParam.output='fourier';% Output type for the analysis
TFParam.channel='EEG';% Specify EEG channels for analysis
TFParam.method='mtmconvol';% Multi-taper method for time-frequency analysis
TFParam.taper='hanning';% Type of tapering applied
TFParam.foi=1:2:30;% Frequencies of interest (1-30 Hz in steps of 2 Hz)
TFParam.t_ftimwin=ones(length(TFParam.foi),1).*0.5;% Time window length (0.5 sec)
TFParam.toi=0:0.5:TimeRange(2)-TimeRange(1);% Time intervals for analysis
%% Main Loop for Subjects
for i=1:length(EDFDir)% Loop through all subjects
    disp(['Subject #' num2str(i)]);
    FileName=[EEGDataPath EDFDir(i).name];
    EEG=pop_biosig(FileName);
    EEG=eeg_checkset(EEG);
    %[ALLEEG,EEG,CURRENTSET,com]=pop_newset(ALLEEG,EEG,0,'gui','off');% Create a new dataset

    %% Create Channel Location
    Elec=readtable([EEGDataPath 'BC-32-X4.txt']);
    EEGChanLoc=FnEEGChanLocCreate(Elec);
    EEG.chanlocs=EEGChanLoc;

    %% PSD Analysis
    Fs=EEG.srate;

    %% Create EEG (32 Chan.) Layout
    OurLayout=FnEEGLayoutCreate(EEGChanLoc);

    %% Read and Add Markers
    FilenameMarkers=[FileName(1:end-4) '.Markers'];
    [EEG,FinalEventName,FinalEventTimes]=FnAddMarkers(EEG,FilenameMarkers);

    %% Epoching Using EEGLAB
    EEG_S14=pop_epoch(EEG,{'S 14'},TimeRange,'newname','BDF file resampled epochs','epochinfo','yes');
    EEG_S16=pop_epoch(EEG,{'S 16'},TimeRange,'newname','BDF file resampled epochs','epochinfo','yes');

    %% Convert EEGLAB Data to FieldTrip Format
    DataS14=eeglab2fieldtrip(EEG_S14,'preprocessing','none');
    DataS16=eeglab2fieldtrip(EEG_S16,'preprocessing','none');

    %% Concatenate Data for Grand-Grand Average TF Analysis
    cfg=[];
    MergedData=ft_appenddata(cfg,DataS14,DataS16);

    %% Apply Time-Frequency Analysis and Baseline Normalization
    TFCat=FnTimeFreqAnalysis(MergedData,TFParam);
    TFCat.dimord='chan_freq_time';
    PowerTemp=FnInducedPower(TFCat.fourierspctrm);
    TFCat.powspctrm=squeeze(nanmean(PowerTemp,1));% Add power data for baseline normalization
    TFCat_BaseNormed=TFCat;

    %% Baseline Normalization and Mean Calculation in Baseline Window
    cfg=[];cfg.baseline=BaselineWin;% Baseline window (in seconds)
    cfg.parameter={'powspctrm'};% Specify the parameter to normalize
    cfg.baselinetype=BaseLineNormtype;% Baseline type ('db', 'relchange', etc.)
    [TFCat_BaseNormed,meanVals]=ft_freqbaseline_EE(cfg,TFCat);
    %% Time-Frequency Analysis on S14 and S16 EEG Data (Conditions: -5dB and 5dB)
    % This loop performs time-frequency analysis for two EEG conditions: S14 (-5dB) and S16 (5dB).
    % For each condition, we compute induced power, apply baseline normalization,
    % and store both condition-specific and average-normalized power spectra.

    Data_Cond={DataS14, DataS16}; % Cell array holding EEG data for conditions S14 and S16
    for CondNum=1:2 % Loop through the two conditions (S14 and S16)

        % Select EEG data for the current condition
        Data2TF=Data_Cond{CondNum};

        % Apply time-frequency analysis
        TF_Cond=FnTimeFreqAnalysis(Data2TF, TFParam); % Perform time-frequency decomposition
        TF_Cond.dimord='chan_freq_time'; % Specify dimensional order: channels, frequencies, time
        TF_Cond_Fourier=TF_Cond;


        % Compute induced power (unnormalized)
        PowerTemp=FnInducedPower(TF_Cond.fourierspctrm); % Extract induced power from Fourier coefficients
        Powspctrm_NoNorm=squeeze(nanmean(PowerTemp, 1)); % Compute mean power across trials without normalization
        TF_Cond.powspctrm=Powspctrm_NoNorm; % Assign unnormalized power spectrum to the TF_Cond structure

        % Configure baseline normalization
        cfg=[];
        cfg.baseline=BaselineWin; % Define the baseline time window (e.g., [0.5 0.8] seconds)
        cfg.parameter={'powspctrm'}; % Specify the parameter for normalization (e.g., power spectrum)
        cfg.baselinetype=BaseLineNormtype; % Define normalization type: 'db', 'relchange', etc.

        % Apply baseline normalization (condition-specific)
        TF_Cond=ft_freqbaseline(cfg, TF_Cond); % Normalize power spectrum using the baseline window
        Powspctrm_CondNorm=TF_Cond.powspctrm; % Extract the condition-specific baseline-normalized power spectrum

        % Compute average baseline normalization
        Powspctrm_AvgNorm=10 * log10(Powspctrm_NoNorm ./ meanVals); % Normalize power spectrum to average baseline values

        % Store results for comparison
        TF_Cond.powspctrm_CondNorm = Powspctrm_CondNorm; % Store condition-specific baseline-normalized power
        TF_Cond.powspctrm_AvgNorm = Powspctrm_AvgNorm; % Store average baseline-normalized power for cross-condition comparisons
        TF_Cond.powspctrm = Powspctrm_NoNorm; % Retain the unnormalized power spectrum as the default representation

        % Calculate coherent power using eigen decomposition of coherence power
        % At this stage, I apply the cross-spectrum, averaged over time within a specific time window.
        % - Time resolution: 0.5-second steps
        % - Frequency resolution: 2 Hz
        % This step focuses on extracting coherence information for specific frequency bands.
        % The analysis is demonstrated specifically for the alpha band (e.g., 8-13 Hz).

        % Reshape Fourier spectrum to dimensions: channels x trials x frequencies x time
        Xw=permute(TF_Cond_Fourier.fourierspctrm,[2,1,3,4]); % Reordering dimensions for further analysis

        % Perform eigen decomposition of the cross-spectra to compute coherence power
        [XYw,Cmat,Ctot,Cvec,Cent,SDiag,SDiagAvg]=FnEigCrossSpectrum(Xw); % Calculate cross-spectra
        SDiagAvg_perm=permute(SDiagAvg,[3,1,2]); % Permute dimensions for further analysis
        CohPow=permute(abs(Cvec),[3,1,2]); % Compute coherence power

        % Loop through each time range and calculate mean power for various normalization types
        %figure,
        for TimeRangeInd=2:size(TimeWindowMatrix,1)
            MeanPower_CondNorm(:,TimeRangeInd,CondNum,i)=FnFindInducedPowerinSelTime(FreqRange,TimeWindowMatrix(TimeRangeInd,:),TF_Cond,TF_Cond.powspctrm_CondNorm);
            MeanPower_AvgNorm(:,TimeRangeInd,CondNum,i)=FnFindInducedPowerinSelTime(FreqRange,TimeWindowMatrix(TimeRangeInd,:),TF_Cond,TF_Cond.powspctrm_AvgNorm);
            MeanPower_NoNorm(:,TimeRangeInd,CondNum,i)=FnFindInducedPowerinSelTime(FreqRange,TimeWindowMatrix(TimeRangeInd,:),TF_Cond,TF_Cond.powspctrm);
            MeanPower_SDiag(:,TimeRangeInd,CondNum,i)=FnFindInducedPowerinSelTime(FreqRange,TimeWindowMatrix(TimeRangeInd,:),TF_Cond,SDiagAvg_perm);
            MeanPower_ChoPow(:,TimeRangeInd,CondNum,i)=FnFindInducedPowerinSelTime(FreqRange,TimeWindowMatrix(TimeRangeInd,:),TF_Cond,CohPow);

            % subplot(4,3,TimeRangeInd)
            % topoplot(MeanPower_CondNorm(:,TimeRangeInd),EEGChanLoc,'electrodes','off','style','both','plotrad',.7,'headrad',.66); % Topographic plot for -5 dB condition
            % colorbar;
        end


        % Method 2: My Ccustom implementation
        % Use a multi-taper FFT approach and segment the data manually for coherence analysis.

    end

end

%% Grand average calculation
AllPower(:,:,:,1)=nanmean(MeanPower_CondNorm,4);
AllPower(:,:,:,2)=nanmean(MeanPower_AvgNorm,4);
AllPower(:,:,:,3)=nanmean(MeanPower_NoNorm,4);
AllPower(:,:,:,4)=nanmean(MeanPower_SDiag,4);
AllPower(:,:,:,5)=nanmean(MeanPower_ChoPow,4);

% Number of conditions, methods, and time windows
NumConditions=2;  % Two conditions (e.g., -5 dB and 5 dB)
NumMethod=5;  % Five different normalization methods
NumWindow=size(TimeWindowMatrix,1);  % Number of time windows (e.g., 12 windows)

% Loop through each condition
for Cond=1:NumConditions
    
    % Extract power for the current condition
    CurrentPow=[];
    CurrentPow=squeeze(AllPower(:,:,Cond,:));

    % Create a new figure for each condition
    figure;
    title(['HC-AvgMast-Cond ' num2str(Cond)]);
    for MethodIdx=1:NumMethod  % Loop through each normalization method
        for WindowIdx=2:NumWindow  % Loop through each time window
            % Determine subplot position for the current method and window
            subplot(NumMethod,NumWindow,(MethodIdx-1)*NumWindow+WindowIdx);

            % Extract the power data for the current method and window
            CurrentPowTemp=[];
            CurrentPowTemp=CurrentPow(:,WindowIdx,MethodIdx);

            % Create topographic plot for the current window and method
            topoplot(CurrentPowTemp,EEGChanLoc,'electrodes','off','style','both','plotrad',.7,'headrad',.66); % Topographic plot for condition
            clim([min(CurrentPowTemp),max(CurrentPowTemp)]);

            % % Check if the brewermap toolbox is available; if not, it will load it
            % ft_hastoolbox('brewermap', 1);  % Ensure that the brewermap toolbox is on the path
            % % Change the colormap for the topographic plots to 'RdBu' and reverse the color order
            % % This colormap is suitable for visualizing EEG microstates with distinct colors.
            % colormap(flipud(brewermap(64, 'RdBu')));  % Set the colormap to RdBu and flip it for better visual contrast
        end
    end
end




