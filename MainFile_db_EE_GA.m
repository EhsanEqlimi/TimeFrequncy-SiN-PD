% This main file performs time-frequency analysis on continnous EEG signals
% Ehsan Eqlimi, @WAVES, UGent,Belgium November 2020
clc;clear;close all;
%% Initialization
warning('off');
currentFolder = pwd;
EEGDataPath=[currentFolder '\Data_PP_Corrected\'];
addpath([currentFolder '\fieldtrip-20200220'])% feildtrip toolbox
addpath([currentFolder '\eeglab2019_1']); %Addpath EEGlab toolbox for filtering, reading loc file (electrode), and plotting PSD topographic maps
addpath([currentFolder '\FASTER']); %Addpath EEGlab toolbox for filtering, reading loc file (electrode), and plotting PSD topographic maps
%Note: BioSig plugin is needed. BioSig is a plugin for eeglab to read EDF
%data.
% Determine where your m-file's folder is.
% folder = fileparts(which('MainFile.m'));
% % Add that folder plus all subfolders to the path.
% addpath(genpath(folder));
Domain='avgall';
EDFDir=dir(fullfile(EEGDataPath,[ '*' Domain '.edf']));%'*corrected.edf'
[ALLEEG, ~, ~, ~] = eeglab;
%% PSD Parameters
Frames=0;
SelFreq=[2,5,10,20,40]; % Selected frequency
Lim=[0 60 -20 70 -10 10]; % Limit
Nchannel=32;
%% Epoching parameters
TimeRange=[0 6]; %Second or [0 4.5]
%% Time frequency parameters
TFParam.pad =[];
TFParam.keeptrials= 'yes';
TFParam.output= 'fourier';
TFParam.channel= 'EEG';
TFParam.method='mtmconvol';
TFParam.taper='hanning';
TFParam.foi= 1:2:30;% analysis 1 to 13 Hz in steps of 2 Hz
TFParam.t_ftimwin=ones(length(TFParam.foi),1).*0.5; %length of time window = 0.5 sec
TFParam.toi= 0:0.05:TimeRange(2)-TimeRange(1);
Method='plv';
%% Main loop for subjects
for i=1:length(EDFDir) % i-->Subject
    disp(['Subject #' num2str(i)]);
    FileName=[EEGDataPath EDFDir(i).name];
    EEG=pop_biosig(FileName);
    EEG=eeg_checkset(EEG);
    [ALLEEG, EEG, CURRENTSET, com]=pop_newset(ALLEEG, EEG, 0,'gui','off'); % And make this a new set
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
    cfg=[]; % create Fieldtrip-like structure
    cfg.datachan=1:size(EEG.data,1); % select scalp electrodes
    % cfg.eyechan = [65 66 67 68]; % select eye channels (useful if running ICA)
    cfg.thresh=[3 3 3 3 3 12]; % see help eegF_FASTER for a description of each number. Lower numbers are more conservative.
    trials2remove14=[];
    [gen_bad_chans,EEG_S14,trials2remove14]=eegF_FASTER_OnlyEpoch(cfg,EEG_S14); % run eegF_FASTER function
    EEG_S14=pop_select(EEG_S14,'notrial',trials2remove14); % remove bad epochs
    trials2remove16=[];
    [gen_bad_chans,EEG_S16,trials2remove16]=eegF_FASTER_OnlyEpoch(cfg,EEG_S16); % run eegF_FASTER function
    EEG_S16=pop_select(EEG_S16,'notrial',trials2remove16); % remove bad epochs
    %% Converting EEGlab data to filedtrip data
    DataS14 = eeglab2fieldtrip(EEG_S14,'preprocessing','none');
    DataS16 = eeglab2fieldtrip(EEG_S16,'preprocessing','none');
    %% Concatenating S14 and S16 for finding Grand-Grand average TF
    cfg = [];
    MergedData =  ft_appenddata(cfg, DataS14, DataS16);
    %% Apply TF on Merged Data and baseline normalization
    TFCat=FnTimeFreqAnalysis(MergedData,TFParam);
    TFCat.dimord='rpt_chan_freq_time';
    cfg=[];  cfg.baseline=[0.5 0.8]; %second
    cfg.parameter={'fourierspctrm'};
    [TFCat_BaseNormed] = ft_freqbaseline(cfg, TFCat); %Baseline normalization
    InducedPowerCat=FnInducedPower(TFCat_BaseNormed.fourierspctrm);%Induced Power
    MeanIndPowerTrialsCat=squeeze(nanmean(InducedPowerCat,1));%Average over trilas
    MeanIndPowerTrialsChannCat(:,:,i)=squeeze(nanmean(MeanIndPowerTrialsCat,1));%Average over channels
    [EvokedPowerCat,MeanEvokedTrialsChannCat(:,:,i)]=FnEvokedPower(TFCat_BaseNormed.fourierspctrm,Method);%Evoked power
    if 0
        %% Time-Frequency analysis
        if 0
            %*********************** Create S 14 (-5 dB SNR) **********************
            DataS14=FnCreateEpochedData(FileName,FinalEventName,FinalEventTimes,TimeRange,Fs,'S 14');
            %****************** Create S16 (+5 dB SNR)**************************
            DataS16=FnCreateEpochedData(FileName,FinalEventName,FinalEventTimes,TimeRange,Fs,'S 16');
        end
        %******************Time-Frequncy for  S14 ****************************
        TFS14=FnTimeFreqAnalysis(DataS14,TFParam);
        TFS14.dimord='rpt_chan_freq_time';
        cfg=[];  cfg.baseline=[0.5 0.8]; %second
        cfg.parameter={'fourierspctrm'};
        [TFS14] = ft_freqbaseline(cfg, TFS14); %Baseline normalization
        %**************************Time-Frequncy for  S16 *********************
        TFS16=FnTimeFreqAnalysis(DataS16,TFParam);
        TFS16.dimord='rpt_chan_freq_time';
        cfg=[];  cfg.baseline=[0.5 0.8]; %second
        cfg.parameter={'fourierspctrm'};
        [TFS16] = ft_freqbaseline(cfg, TFS16); %Baseline normalization
        %*********************** Induced Power for S14 ************************
        PowerS14=FnInducedPower(TFS14.fourierspctrm);
        MeanPowerTrialsS14=squeeze(nanmean(PowerS14,1));
        [SelLab,SelInd,~]=intersect(TFS14.label,{'F3','Fz','F4'},'stable'); % select electrodes for induced
        MeanPowerS14=squeeze(nanmean(MeanPowerTrialsS14(SelInd,:,:),1));
        %*********************** Induced Power for S16 ***********************
        PowerS16=FnInducedPower(TFS16.fourierspctrm);
        MeanPowerTrialsS16=squeeze(nanmean(PowerS16,1));
        MeanPowerS16=squeeze(nanmean(MeanPowerTrialsS16(SelInd,:,:),1));
        %************************************ Plot S14 induced power ************
        figure,surf(TFParam.toi,TFParam.foi,MeanPowerS14);
        caxis([min(MeanPowerS14(:)) max(MeanPowerS14(:))])
        shading interp;
        xlabel('Time (sec)');ylabel('Frequency (Hz)');zlabel('Power spectrum');
        view(0,-90)
        h=colorbar;
        colormap('jet');
        ylabel(h, 'Power averaged across all channels for the whole trial duration')
        title(['-5 dB SNR, S14-Subject-' num2str(i)])
        %************************************ Plot S16 induced power ************
        figure,surf(TFParam.toi,TFParam.foi,MeanPowerS16);
        caxis([min(MeanPowerS14(:)) max(MeanPowerS14(:))])
        shading interp;
        xlabel('Time (sec)');ylabel('Frequency (Hz)');zlabel('Power spectrum');
        view(0,-90)
        h=colorbar;
        colormap('jet');
        ylabel(h, 'Power averaged across all channels for the whole trial duration')
        title(['+5 dB SNR, S16-Subject-' num2str(i)])
        %% Evoked Power
        %**************************** Evoked Power S14 ************************
        [EvokedPowerS14,MeanEvokedPowerS14]=FnEvokedPower(TFS14.fourierspctrm,Method);
        %**************************** Evoked Power S16 ************************
        [EvokedPowerS16,MeanEvokedPowerS16]=FnEvokedPower(TFS16.fourierspctrm,Method);
        %*******************************Plot S14 evoked power******************
        figure,surf(TFParam.toi,TFParam.foi,MeanEvokedPowerS14);
        caxis([min(min(MeanEvokedPowerS14(:)),min(MeanEvokedPowerS16(:))) max(max(MeanEvokedPowerS14(:)),max(MeanEvokedPowerS14(:)))])
        shading interp;
        xlabel('Time (sec)');ylabel('Frequency (Hz)');zlabel('Power spectrum');
        view(0,-90)
        h=colorbar;
        colormap('jet');
        ylabel(h, [Method '-Evoked Power averaged across all channels for the whole trial duration']);
        title(['-5 dB SNR, S14-Subject-' num2str(i)])
        %*******************************Plot S16 evoked power******************
        figure,surf(TFParam.toi,TFParam.foi,MeanEvokedPowerS16);
        caxis([min(min(MeanEvokedPowerS14(:)),min(MeanEvokedPowerS16(:))) max(max(MeanEvokedPowerS14(:)),max(MeanEvokedPowerS14(:)))])
        shading interp;
        xlabel('Time (sec)');ylabel('Frequency (Hz)');zlabel('Power spectrum');
        view(0,-90)
        h=colorbar;
        colormap('jet');
        ylabel(h, [Method '-Evoked Power averaged across all channels for the whole trial duration']);
        title(['+5 dB SNR, S16-Subject-' num2str(i)])
        % ***************************** pre-stimulus induced alpha-S 14*************
        FreqRange=[8 12]; %Hz
        TimeRange_PreSent=[1.5 2]; %Second
        MeanInducedPowerS14_PreSent=FnFindInducedPowerinSelTime(FreqRange,TimeRange_PreSent,TFS14,MeanPowerTrialsS14);
        
        TimeRange_DuringSent=[2.5 4.5]; %Second
        MeanInducedPowerS14_DuringSent=FnFindInducedPowerinSelTime(FreqRange,TimeRange_DuringSent,TFS14,MeanPowerTrialsS14);
        
        TimeRange_PostSent=[4.5 6]; %Second
        MeanInducedPowerS14_PostSent=FnFindInducedPowerinSelTime(FreqRange,TimeRange_PostSent,TFS14,MeanPowerTrialsS14);
        % ***************************** Plot pre-stimulus induced alpha- S 14************
        Type='Induced Alpha Power';
        FnTopoPlotPower(MeanInducedPowerS14_PreSent,EEGChanLoc,0,max(MeanInducedPowerS14_PreSent),Type);
        title(['-5 dB SNR, S14-Subject-' num2str(i)  '-[' num2str(TimeRange_PreSent(1)) '-'...
            num2str(TimeRange_PreSent(2)) '] sec-[' num2str(FreqRange(1)) '-' num2str(FreqRange(2)) '] Hz']);
        
        FnTopoPlotPower(MeanInducedPowerS14_DuringSent,EEGChanLoc,0,max(MeanInducedPowerS14_DuringSent),Type);
        title(['-5 dB SNR, S14-Subject-' num2str(i)  '-[' num2str(TimeRange_DuringSent(1)) '-'...
            num2str(TimeRange_DuringSent(2)) '] sec-[' num2str(FreqRange(1)) '-' num2str(FreqRange(2)) '] Hz']);
        
        FnTopoPlotPower(MeanInducedPowerS14_PostSent,EEGChanLoc,0,max(MeanInducedPowerS14_PostSent),Type);
        title(['-5 dB SNR, S14-Subject-' num2str(i)  '-[' num2str(TimeRange_PostSent(1)) '-'...
            num2str(TimeRange_PostSent(2)) '] sec-[' num2str(FreqRange(1)) '-' num2str(FreqRange(2)) '] Hz']);
        % ***************************** pre-stimulus induced alpha-S 16*************
        FreqRange=[8 12]; %Hz
        TimeRange_PreSent=[1.5 2]; %Second
        MeanInducedPowerS16_PreSent=FnFindInducedPowerinSelTime(FreqRange,TimeRange_PreSent,TFS16,MeanPowerTrialsS16);
        
        TimeRange_DuringSent=[2.5 4.5]; %Second
        MeanInducedPowerS16_DuringSent=FnFindInducedPowerinSelTime(FreqRange,TimeRange_DuringSent,TFS16,MeanPowerTrialsS16);
        
        TimeRange_PostSent=[4.5 6]; %Second
        MeanInducedPowerS16_PostSent=FnFindInducedPowerinSelTime(FreqRange,TimeRange_PostSent,TFS16,MeanPowerTrialsS16);
        % ***************************** Plot pre-stimulus induced alpha- S 14************
        FnTopoPlotPower(MeanInducedPowerS16_PreSent,EEGChanLoc,0,max(MeanInducedPowerS14_PreSent),Type);
        title(['+5 dB SNR, S16-Subject-' num2str(i)  '-[' num2str(TimeRange_PreSent(1)) '-'...
            num2str(TimeRange_PreSent(2)) '] sec-[' num2str(FreqRange(1)) '-' num2str(FreqRange(2)) '] Hz']);
        
        FnTopoPlotPower(MeanInducedPowerS16_DuringSent,EEGChanLoc,0,max(MeanInducedPowerS14_DuringSent),Type);
        title(['+5 dB SNR, S16-Subject-' num2str(i)  '-[' num2str(TimeRange_DuringSent(1)) '-'...
            num2str(TimeRange_DuringSent(2)) '] sec-[' num2str(FreqRange(1)) '-' num2str(FreqRange(2)) '] Hz']);
        
        FnTopoPlotPower(MeanInducedPowerS16_PostSent,EEGChanLoc,0,max(MeanInducedPowerS14_PostSent),Type);
        title(['+5 dB SNR, S16-Subject-' num2str(i)  '-[' num2str(TimeRange_PostSent(1)) '-'...
            num2str(TimeRange_PostSent(2)) '] sec-[' num2str(FreqRange(1)) '-' num2str(FreqRange(2)) '] Hz']);
        
        %************************************Multi Plot Power************************
        %     cfg = [];
        %     % cfg.baseline     = [-0.5 -0.1];
        %     cfg.baselinetype = 'absolute';
        %     cfg.zlim         = [min(MeanPowerS14(:)) max(MeanPowerS14(:))];
        %     cfg.showlabels   = 'yes';
        %     cfg.layout       = OurLayout;
        %     cfg.colorbar     = 'yes';
        %     figure
        %     ft_multiplotTFR(cfg, MeanPowerS14)
        %     % colormap('jet');
        %     title(['-5 dB SNR, S14-Subject-' num2str(i) '-Induced Power'])
        
        %************************************Multi Plot PLV************************
        %     cfg = [];
        %     % cfg.baseline     = [-0.5 -0.1];
        %     cfg.baselinetype = 'absolute';
        %     cfg.zlim         = [min(MeanPLV_S14(:)) max(MeanPLV_S14(:))];
        %     cfg.showlabels   = 'yes';
        %     cfg.layout       = OurLayout;
        %     cfg.colorbar     = 'yes';
        %     figure
        %     Freq_S14_PLV=Freq_S14;
        %     Freq_S14_PLV.powspctrm=permute(PlVTap,[2 1 3]);
        %     ft_multiplotTFR(cfg, Freq_S14_PLV)
        %     % colormap('jet');
        %     title(['-5 dB SNR, S14-Subject-' num2str(i) '-PLV'])
        %% Post-MTB evoked power
        % ***************************** post-MTB evoked  LF power S14****************
        FreqRange=[2 8]; %Hz
        TimeRange_PostMTB=[1 1.5]; %Second
        MeanEvokedPowerS14_PostMTB=FnFindEvokedPowerinSelTime(FreqRange,TimeRange_PostMTB,TFS14,EvokedPowerS14);
        
        TimeRange_PostSent=[2 2.5]; %Second
        MeanEvokedPowerS14_PostSent=FnFindEvokedPowerinSelTime(FreqRange,TimeRange_PostSent,TFS14,EvokedPowerS14);
        % ***************************** Post-MTB evoked  LF power S16*************
        FreqRange=[2 8]; %Hz
        TimeRange_PostMTB=[1 1.5]; %Second
        MeanEvokedPowerS16_PostMTB=FnFindEvokedPowerinSelTime(FreqRange,TimeRange_PostMTB,TFS16,EvokedPowerS16);
        
        TimeRange_PostSent=[2 2.5]; %Second
        MeanEvokedPowerS16_PostSent=FnFindEvokedPowerinSelTime(FreqRange,TimeRange_PostSent,TFS16,EvokedPowerS16);
        % ***************************** Plot evoked alpha- S 14************
        Type='Evoked Low-Freq Power';
        FnTopoPlotPower(MeanEvokedPowerS14_PostMTB,EEGChanLoc,min(min(MeanEvokedPowerS14_PostMTB),min(MeanEvokedPowerS16_PostMTB)),max(max(MeanEvokedPowerS14_PostMTB),max(MeanEvokedPowerS16_PostMTB)),Type);
        title(['-5 dB SNR, S14-Subject-' num2str(i)  '-[' num2str(TimeRange_PostMTB(1)) '-'...
            num2str(TimeRange_PostMTB(2)) '] sec-[' num2str(FreqRange(1)) '-' num2str(FreqRange(2)) '] Hz']);
        
        FnTopoPlotPower(MeanEvokedPowerS14_PostSent,EEGChanLoc,min(min(MeanEvokedPowerS14_PostSent),min(MeanEvokedPowerS16_PostSent)),max(max(MeanEvokedPowerS14_PostSent),max(MeanEvokedPowerS16_PostSent)),Type);
        title(['-5 dB SNR, S14-Subject-' num2str(i)  '-[' num2str(TimeRange_PostSent(1)) '-'...
            num2str(TimeRange_PostSent(2)) '] sec-[' num2str(FreqRange(1)) '-' num2str(FreqRange(2)) '] Hz']);
        % ***************************** Plot evoked alpha- S 16************
        FnTopoPlotPower(MeanEvokedPowerS16_PostMTB,EEGChanLoc,min(min(MeanEvokedPowerS14_PostMTB),min(MeanEvokedPowerS16_PostMTB)),max(max(MeanEvokedPowerS14_PostMTB),max(MeanEvokedPowerS16_PostMTB)),Type);
        title(['+5 dB SNR, S16-Subject-' num2str(i)  '-[' num2str(TimeRange_PostMTB(1)) '-'...
            num2str(TimeRange_PostMTB(2)) '] sec-[' num2str(FreqRange(1)) '-' num2str(FreqRange(2)) '] Hz']);
        
        FnTopoPlotPower(MeanEvokedPowerS16_PostSent,EEGChanLoc,min(min(MeanEvokedPowerS14_PostSent),min(MeanEvokedPowerS16_PostSent)),max(max(MeanEvokedPowerS14_PostSent),max(MeanEvokedPowerS16_PostSent)),Type);
        title(['+5 dB SNR, S16-Subject-' num2str(i)  '-[' num2str(TimeRange_PostSent(1)) '-'...
            num2str(TimeRange_PostSent(2)) '] sec-[' num2str(FreqRange(1)) '-' num2str(FreqRange(2)) '] Hz']);
        %% Find Fz index
        [SelLab,SelInd,~]=intersect(TFS14.label,{'F3','Fz','F4'},'stable');
        
        AllPreStimAlphaPower(i,1)=nanmean(MeanInducedPowerS14_PreSent(SelInd));
        AllPreStimAlphaPower(i,2)=nanmean(MeanInducedPowerS16_PreSent(SelInd));
        % Find central channel indices
        [CentrealElec,CentrealElec1,CentrealEle2]=intersect(TFS14.label,{'FC1','FC2','FCz','C3','Cz','C4'},'stable');
        AllPostMTBLFPLV(i,1)=nanmean(MeanEvokedPowerS14_PostMTB(CentrealElec1)); %% also store all electrodes for post-processing
        AllPostMTBLFPLV(i,2)=nanmean(MeanEvokedPowerS16_PostMTB(CentrealElec1));
        %% Writing
        AllPowers=[];
        AllPowers=[MeanInducedPowerS14_PreSent MeanInducedPowerS14_DuringSent MeanInducedPowerS14_PostSent ...
            MeanInducedPowerS16_PreSent MeanInducedPowerS16_DuringSent MeanInducedPowerS16_PostSent...
            MeanEvokedPowerS14_PostMTB MeanEvokedPowerS14_PostSent...
            MeanEvokedPowerS16_PostMTB MeanEvokedPowerS16_PostSent];
        
        VarNames={'Channles','IndS14PreSent','IndS14DurSent','IndS14PostSent',...
            'IndS16PreSent','IndS16DurSent','IndS16PostSent'...
            'EvS14PostMTB','EvS14PostSent',...
            'EvS16PostMTB','EvS16PostSent'
            };%Variable names
        Table = table(TFS14.label,MeanInducedPowerS14_PreSent,MeanInducedPowerS14_DuringSent,MeanInducedPowerS14_PostSent, ...
            MeanInducedPowerS16_PreSent,MeanInducedPowerS16_DuringSent,MeanInducedPowerS16_PostSent,...
            MeanEvokedPowerS14_PostMTB,MeanEvokedPowerS14_PostSent,...
            MeanEvokedPowerS16_PostMTB,MeanEvokedPowerS16_PostSent);%put in a table
        Table.Properties.VariableNames=VarNames;%assign the variable names
        writetable(Table,[EDFDir(i).name(1:end-4) '_TimeFreqPower' '.txt'],'Delimiter',' '); %write
        %% Keepinp the data for Grand-average time-frequaency represneation across trial, subjects, and electrodes
        InducedTFS16(:,:,i)=MeanPowerS16;
        InducedTFS14(:,:,i)=MeanPowerS14;
        EvokedTFS16(:,:,i)=MeanEvokedPowerS16;
        EvokedTFS14(:,:,i)=MeanEvokedPowerS14;
    end
end
%%  Grand-average time-frequaency represneation across trial, conditions, subjects, and electrodes
GrandAvgInducedTF=nanmean(MeanIndPowerTrialsChannCat,3);
GrandAvgEvokedTF=nanmean(MeanEvokedTrialsChannCat,3);

figure,surf(TFParam.toi,TFParam.foi,GrandAvgInducedTF);
% caxis([min(min(GrandAvgInducedTF(:)),min(GrandAvgInducedTF(:))) max(max(GrandAvgInducedTF(:)),max(GrandAvgInducedTF(:)))])
shading interp;
xlabel('Time (sec)');ylabel('Frequency (Hz)');zlabel('Power spectrum');
view(0,-90)
h=colorbar;
colormap('jet');
ylabel(h, ['Induced Grand-average TF across trial, subjects, conditions and electrodes']);


figure,surf(TFParam.toi,TFParam.foi,GrandAvgEvokedTF);
caxis([min(min(GrandAvgEvokedTF(:)),min(GrandAvgEvokedTF(:))) max(max(GrandAvgEvokedTF(:)),max(GrandAvgEvokedTF(:)))])
shading interp;
xlabel('Time (sec)');ylabel('Frequency (Hz)');zlabel('Power spectrum');
view(0,-90)
h=colorbar;
colormap('jet');
ylabel(h, ['Evoked Grand-average TF across trial, subjects, conditions and electrodes']);
if 0
    %% plotting overall results
    
    XVals=[1:size(AllPreStimAlphaPower,1)];
    figure,bar(XVals,AllPreStimAlphaPower)
    xlabel('Subject No.');
    xticks([1 2 3]);
    xticklabels({'HC1', 'HC2', 'PD1'});
    ylabel('Pre-sentence Alpha Induced Power at Fz');
    legend({'-5 dB','+5 dB'},'location','Best');
    
    figure,bar(XVals,AllPostMTBLFPLV(:,1))
    bar(XVals,AllPostMTBLFPLV);
    xlabel('Subject No.');
    xticks([1 2 3]);
    xticklabels({'HC1', 'HC2', 'PD1'});
    ylabel('Post-MTB LF Evoked Power (PLV) over Centeral Channels');
    legend({'-5 dB','+5 dB'});
end
%% output to txt --> excell
%% time windows : evoked [1 1.5] and [2 2.5]; indiced alpha [1.5 2] [2.5 4.5] [4.5 6]

