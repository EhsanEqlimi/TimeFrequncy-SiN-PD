% This main file performs time-frequency analysis on continnous EEG signals
% Ehsan Eqlimi, @WAVES, UGent,Belgium November 2020
clc;
% clear;
close all;
%% Initialization
warning('off');
currentFolder = pwd;


% EEGDataPath=[currentFolder '\Data_PP_Corrected\'];
% addpath([currentFolder '\fieldtrip-20200220'])% feildtrip toolbox
% addpath([currentFolder '\eeglab2019_1']); %Addpath EEGlab toolbox for filtering, reading loc file (electrode), and plotting PSD topographic maps
% addpath([currentFolder '\FASTER']); %Addpath EEGlab toolbox for filtering, reading loc file (electrode), and plotting PSD topographic maps
%Note: BioSig plugin is needed. BioSig is a plugin for eeglab to read EDF
%data.
% Determine where your m-file's folder is.
% folder = fileparts(which('MainFile.m'));
% % Add that folder plus all subfolders to the path.
% addpath(genpath(folder));
EEGDataPath='E:\Ehsan\Data_PP_Corrected\';
%% Add Required Toolboxes to Path
% Add the FieldTrip toolbox for EEG analysis
addpath('E:\Ehsan\RippleServerFiles\Toolbox\fieldtrip-20220729'); % Add FieldTrip to path
ft_defaults; % Initialize FieldTrip



Domain='PD*avgall';%'avgmast'; %'avgall';
EDFDir=dir(fullfile(EEGDataPath,[ '*' Domain '.edf']));%'*corrected.edf'
% [ALLEEG, ~, ~, ~] = eeglab;
% eeglab
%% PSD Parameters
Frames=0;
SelFreq=[2,5,10,20,40]; % Selected frequency
Lim=[0 60 -20 70 -10 10]; % Limit
Nchannel=32;
%% Epoching parameters
TimeRange=[0 6]; %Second or [0 4.5]
%% Time frequency parameters
% TFParam.pad =[];
% TFParam.keeptrials= 'yes';
% TFParam.keeptapers='yes';
% TFParam.output= 'fourier';%'fourier';
% TFParam.channel= 'EEG';
% TFParam.method='mtmconvol';%'mtmfft';%'mtmconvol';
% TFParam.taper='dpss';%'dpss';%'hanning';
% TFParam.foi= 1:1:30;% analysis 1 to 13 Hz in steps of 2 Hz
% TFParam.tapsmofrq=1;
% TFParam.t_ftimwin=ones(length(TFParam.foi),1).*1; %length of time window = 0.5 sec
% TFParam.toi= 0:0.05:TimeRange(2)-TimeRange(1);
% Method='itpc';



TFParam.pad =[];
TFParam.keeptrials= 'yes';
TFParam.output= 'fourier';
TFParam.channel= 'EEG';
TFParam.method='mtmconvol';
TFParam.taper='hanning';
TFParam.foi= 1:2:30;% analysis 1 to 13 Hz in steps of 2 Hz
TFParam.t_ftimwin=ones(length(TFParam.foi),1).*0.5; %length of time window = 0.5 sec
TFParam.toi= 0:0.05:TimeRange(2)-TimeRange(1);
Method='itpc';
SegmentedTFCat=[];
%% Main loop for subjects
for i=[1 3:length(EDFDir)] % i-->Subject
    disp(['Subject #' num2str(i)]);
    FileName=[EEGDataPath EDFDir(i).name];
    EEG=pop_biosig(FileName);
    EEG=eeg_checkset(EEG);
    % [ALLEEG, EEG, CURRENTSET, com]=pop_newset(ALLEEG, EEG, 0,'gui','off'); % And make this a new set
    %%  Create Channel location
    Elec=readtable([EEGDataPath 'BC-32-X4.txt']);
    EEGChanLoc=FnEEGChanLocCreate(Elec);
    EEG.chanlocs=EEGChanLoc;
    %% PSD
    Fs=EEG.srate;
    if 0
        Title=['PSD-' EDFDir(i).name];
        figure,
        [Spectra,Freqs,Speccomp,Ccontrib,Specstd] = ...
            spectopo(EEG.data, Frames, Fs, 'freq',SelFreq,...
            'chanlocs',EEGChanLoc,'limits',Lim,'title', Title,'electrodes','labels');
    end
    %% If you want to reject bad and noisy EEg channels, uncomment below:
    % [EEG,indelec] = pop_rejchan(EEG, 'elec',[1:32] ,'threshold',5,'norm','on','measure','kurt');
    %% Create EEG (32 Chan.) Layout (it is needed for multi-topoplot)
    OurLayout=FnEEGLayoutCreate(EEGChanLoc);
    %% Read and add markers
    %%%% Read marker file
    FilenameMarkers=[FileName(1:end-4) '.Markers'];
    [EEG,FinalEventName,FinalEventTimes]=FnAddMarkers(EEG,FilenameMarkers);
    %% Epoching using EEGLab (optional)
    % EEG_S14 = pop_epoch( EEG, {'S 14'} , [0 6], 'newname', 'BDF file resampled epochs', 'epochinfo', 'yes');
    % EEG_S16 = pop_epoch( EEG, {'S 16'} , [0 6], 'newname', 'BDF file resampled epochs', 'epochinfo', 'yes');
    %% Write channel added EDF file (optional)
    % pop_writeeeg(EEG, [FileName(1:end-4) '_ChannelAdded.edf'], 'TYPE','EDF');
    %% Remove bad epochs (0-4.5 s) using EEGlab
    EEG_S14 = pop_epoch( EEG, {'S 14'} , [0 6], 'newname', 'BDF file resampled epochs', 'epochinfo', 'yes');
    %     EEG_S14 = pop_eegthresh(EEG_S14,1,[1:size(EEG.data,1)] ,-200,200,0,4.5,0,1);%Artifact rejection
    
    EEG_S16 = pop_epoch( EEG, {'S 16'} , [0 6], 'newname', 'BDF file resampled epochs', 'epochinfo', 'yes');
    %     EEG_S16= pop_eegthresh(EEG_S16,1,[1:size(EEG.data,1)] ,-200,200,0,4.5,0,1);%Artifact rejection
    %% Remove bad epochs using FASTER
    % FASTER: Fully Automated Statistical Thresholding for EEG artifact Rejection (Nolan et al., 2010, J Neurosci Methods)
    % cfg=[]; % create Fieldtrip-like structure
    % cfg.datachan=1:size(EEG.data,1); % select scalp electrodes
    % % cfg.eyechan = [65 66 67 68]; % select eye channels (useful if running ICA)
    % cfg.thresh=[3 3 3 3 3 12]; % see help eegF_FASTER for a description of each number. Lower numbers are more conservative.
    % trials2remove14=[];
    % [gen_bad_chans,EEG_S14,trials2remove14]=eegF_FASTER_OnlyEpoch(cfg,EEG_S14); % run eegF_FASTER function
    % EEG_S14=pop_select(EEG_S14,'notrial',trials2remove14); % remove bad epochs
    % trials2remove16=[];
    % [gen_bad_chans,EEG_S16,trials2remove16]=eegF_FASTER_OnlyEpoch(cfg,EEG_S16); % run eegF_FASTER function
    % EEG_S16=pop_select(EEG_S16,'notrial',trials2remove16); % remove bad epochs
    %% Converting EEGlab data to filedtrip data
    DataS14 = eeglab2fieldtrip(EEG_S14,'preprocessing','none');
    DataS16 = eeglab2fieldtrip(EEG_S16,'preprocessing','none');
    NumTr1(i)=length(DataS14.trial);
    NumTr2(i)=length(DataS16.trial);

    %% Concatenating S14 and S16 for finding Grand-Grand average TF
    cfg = [];
    MergedData=ft_appenddata(cfg, DataS14, DataS16);
    %% Apply TF on Merged Data and baseline normalization
    TFCat=FnTimeFreqAnalysis(MergedData,TFParam);
    TFCat.dimord='rpt_chan_freq_time';
    [NTrials,NChans, NFreqs,NTimes]=size(TFCat.fourierspctrm);
    [~, AlphaIdx_Low]=min(abs(TFCat.freq-8));[~, AlphaIdx_Up]=min(abs(TFCat.freq-13));
    AlphaFreqs=TFCat.freq(AlphaIdx_Low:AlphaIdx_Up);
    for tt=0:6
        [~, TimeIdx(tt+1)]=min(abs(TFCat.time-tt));
    end
    for tt=1:length(TimeIdx)-1
        SegmentedTFCat(i,tt,:,:,:,:)=TFCat.fourierspctrm(:,:,:,TimeIdx(tt):TimeIdx(tt+1));
    end
    % nanmean()
    % 
    % topoplot(SegmentedTFCat{1,1},EEGChanLoc,'electrodes','off','style','both','plotrad',.7,'headrad',.66); % Topographic plot for -5 dB condition
    % colorbar;
end
%% Grand avergae
% Identify slices where all elements are zero along the first dimension
SegmentedTFCat2 = SegmentedTFCat(~all(SegmentedTFCat == 0, [2,3,4,5,6]), :, :, :, :, :);
TF_Minus5dB=SegmentedTFCat2(:,:,1:50,:,AlphaIdx_Low:AlphaIdx_Up,:);
TF_Plus5dB=SegmentedTFCat2(:,:,51:100,:,AlphaIdx_Low:AlphaIdx_Up,:);

Mag_TF_Minus5dB=abs(TF_Minus5dB).^2;
Mag_TF_Plus5dB=abs(TF_Plus5dB).^2;




TF_Minus5dB_GA=squeeze(nanmean(Mag_TF_Minus5dB,[ 1 3 5 6 ]));
TF_Plus5dB_GA=squeeze(nanmean(Mag_TF_Plus5dB,[ 1 3 5 6 ]));

for SegIdx=1:size(TF_Minus5dB_GA,1)
figure,topoplot(Mag_TF_Minus5dB(SegIdx,:),EEGChanLoc,'electrodes','off','style','both','plotrad',.7,'headrad',.66); % Topographic plot for -5 dB condition
end

