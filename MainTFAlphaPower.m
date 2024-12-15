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
BaselineWin=[0.3 0.6];% Baseline window in seconds

% Time ranges for induced power analysis
TimeRange_PreSent_Ind=[1.5 2.0];% Previous selection: [1.5 2]
TimeRange_DuringSent_Ind=[2.5 3.5];% Previous selection: [2.5 4.5]
TimeRange_PostSent_Ind=[4.5 5.5];% Induced power analysis window
BaseLineNormtype='db';% Baseline normalization type (e.g., 'db')
%% Time-Frequency Analysis Parameters
TFParam.pad=[];% Zero-padding (if any)
TFParam.keeptrials='yes';% Keep individual trials during analysis
TFParam.output='fourier';% Output type for the analysis
TFParam.channel='EEG';% Specify EEG channels for analysis
TFParam.method='mtmconvol';% Multi-taper method for time-frequency analysis
TFParam.taper='hanning';% Type of tapering applied
TFParam.foi=1:2:30;% Frequencies of interest (1-30 Hz in steps of 2 Hz)
TFParam.t_ftimwin=ones(length(TFParam.foi),1).*0.5;% Time window length (0.5 sec)
TFParam.toi=0:0.05:TimeRange(2)-TimeRange(1);% Time intervals for analysis
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
    
  
    
end
