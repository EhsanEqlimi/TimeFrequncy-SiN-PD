% This main file performs time-frequency analysis on continnous EEG signals
% Ehsan Eqlimi, @WAVES, UGent,Belgium November 2020
%Last Update: 12/03/2021 by EE
%Last Update: 14/02/2021 by EE
clc;
% clear;
% close all;
%% Initialization
warning('off');
currentFolder = pwd;
DoFaster=0;
% EEGDataPath=[currentFolder '\Data_PP_Corrected\'];
% addpath([currentFolder '\fieldtrip-20200220'])% feildtrip toolbox
% addpath([currentFolder '\eeglab2019_1']); %Addpath EEGlab toolbox for filtering, reading loc file (electrode), and plotting PSD topographic maps
% % addpath('D:\eeglab_current\eeglab2019_1');
% addpath([currentFolder '\FASTER']); %Addpath EEGlab toolbox for filtering, reading loc file (electrode), and plotting PSD topographic maps
% %Note: BioSig plugin is needed. BioSig is a plugin for eeglab to read EDF
%data.
% Determine where your m-file's folder is.
% folder = fileparts(which('MainFile.m'));
% % Add that folder plus all subfolders to the path.
% addpath(genpath(folder));
EEGDataPath='E:\Ehsan\Data_PP_Corrected\';
addpath('E:\Ehsan\GitHub\ComplexPCA-PD');
%% Add Required Toolboxes to Path
% Add the FieldTrip toolbox for EEG analysis
addpath('E:\Ehsan\RippleServerFiles\Toolbox\fieldtrip-20220729'); % Add FieldTrip to path
ft_defaults; % Initialize FieldTrip
% addpath(genpath('E:\Ehsan\GitHub\Toolbox\eeglab2019_1'));
Domain= 'HC*avgall';%'HC*avgmast'; %'avgall';
EDFDir=dir(fullfile(EEGDataPath,[ '*' Domain '.edf']));%'*corrected.edf'
% [ALLEEG, ~, ~, ~] = eeglab;
CatTable=[];
%% PSD Parameters
Frames=0;
SelFreq=[2,5,10,20,40]; % Selected frequency
Lim=[0 60 -20 70 -10 10]; % Limit
Nchannel=32;
%% Epoching and windowing parameters
TimeRange=[0 6]; %Second or [0 4.5]
BaselineWin=[0.3 0.6];
%%%%%Induced parameters:
TimeRange_PreSent_Ind= [1.5 2.0]; %previuous selection: [1.5 2]; %Second
TimeRange_DuringSent_Ind=[2.5 3.5];  %previuous selection: [2.5 4.5]; %Second
TimeRange_PostSent_Ind=[4.5 5.5];%Second
%%%%%%%%Evoked Parameter:
TimeRange_PostMTB_Evk=[1 1.4];% Previuous selction: [1 1.5]; %Second
TimeRange_PreSent_Evk=[1.5 2]; % Presentence antipating alpha power?
TimeRange_PostSent_Evk=[2 2.4];%Second
BaseLineNormtype='db';
Norm=0;
SelChanSwitch=0; %only desired channels;
SelChannLabels={'F4','Fz','F3','C4','Cz','C3','FCz','FC5','FC6','FC2','FC1'};% Induced: 'F4','Fz','F3','C4','Cz','C3','P4','Pz','P3' %Evoked: 'F4','Fz','F3','C4','Cz','C3','FCz','FC5','FC6','FC2','FC1'
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
Method='itpcz';%itpc
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
    EEG_S14 = pop_epoch( EEG, {'S 14'} , TimeRange, 'newname', 'BDF file resampled epochs', 'epochinfo', 'yes');
    %     EEG_S14 = pop_eegthresh(EEG_S14,1,[1:size(EEG.data,1)] ,-200,200,0,4.5,0,1);%Artifact rejection
    
    EEG_S16 = pop_epoch( EEG, {'S 16'} , TimeRange, 'newname', 'BDF file resampled epochs', 'epochinfo', 'yes');
    %     EEG_S16= pop_eegthresh(EEG_S16,1,[1:size(EEG.data,1)] ,-200,200,0,4.5,0,1);%Artifact rejection
    if DoFaster
    %% Remove bad epochs using FASTER
    % FASTER: Fully Automated Statistical Thresholding for EEG artifact Rejection (Nolan et al., 2010, J Neurosci Methods)
    cfg=[]; % create Fieldtrip-like structure
    cfg.datachan=1:size(EEG.data,1); % select scalp electrodes
    % cfg.eyechan = [65 66 67 68]; % select eye channels (useful if running ICA)
    cfg.thresh=[3 3 3 3 3 12]; % see help eegF_FASTER for a description of each number. Lower numbers are more conservative.
    trials2remove14=[];
    [~,EEG_S14,trials2remove14]=eegF_FASTER_OnlyEpoch(cfg,EEG_S14); % run eegF_FASTER function
    EEG_S14=pop_select(EEG_S14,'notrial',trials2remove14); % remove bad epochs
    trials2remove16=[];
    [~,EEG_S16,trials2remove16]=eegF_FASTER_OnlyEpoch(cfg,EEG_S16); % run eegF_FASTER function
    EEG_S16=pop_select(EEG_S16,'notrial',trials2remove16); % remove bad epochs
    end
    %% Converting EEGlab data to filedtrip data
    DataS14 = eeglab2fieldtrip(EEG_S14,'preprocessing','none');
    DataS16 = eeglab2fieldtrip(EEG_S16,'preprocessing','none');
    if 1 %use 1 for Grand averageing
        %% Concatenating S14 and S16 for finding Grand-Grand average TF
        cfg = [];
        MergedData=ft_appenddata(cfg, DataS14, DataS16);
        %% Apply TF on Merged Data and baseline normalization
        TFCat=FnTimeFreqAnalysis(MergedData,TFParam);
        TFCat.dimord='chan_freq_time';
        PowerTemp=FnInducedPower(TFCat.fourierspctrm);
        TFCat.powspctrm=squeeze(nanmean(PowerTemp,1));% I added power data for baseline norm.
        TFCat_BaseNormed=TFCat;
        if Norm
            cfg=[];  cfg.baseline=BaselineWin; %previuous selection: [0.5 0.8]; %second
            cfg.parameter={'powspctrm'}; %fourierspctrm
            cfg.baselinetype=BaseLineNormtype; %or relchange %or dB
            
            [TFCat_BaseNormed,meanVals] = ft_freqbaseline_EE(cfg, TFCat); %Baseline normalization, %meanVals--> baseline
        end
        %%InducedPowerCat=FnInducedPower(TFCat_BaseNormed.fourierspctrm);%Induced Power
        MeanIndPowerTrialsCat=  TFCat_BaseNormed.powspctrm;%Average over trilas
        MeanIndPowerTrialsChannCat(:,:,i)=squeeze(nanmean(MeanIndPowerTrialsCat,1));%Average over channels
        [EvokedPowerCat,MeanEvokedTrialsChannCat(:,:,i)]=FnEvokedPower(TFCat.fourierspctrm,Method);%Evoked power (for evoked power (plv), I did'nt apply baseline normalization, it does not make sense
        %% Find evoked and induced power for merged data
        FreqRange=[8 12]; %Hz
        %%TimeRange_PreSent= [1.5 2.0]; %previuous selection: [1.5 2]; %Second
        MeanInducedPowerCat_PreSent(:,i)=FnFindInducedPowerinSelTime(FreqRange,TimeRange_PreSent_Ind,TFCat_BaseNormed,MeanIndPowerTrialsCat);
        
        %%TimeRange_DuringSent=[2.5 3.5];  %previuous selection: [2.5 4.5]; %Second
        MeanInducedPowerCat_DuringSent(:,i)=FnFindInducedPowerinSelTime(FreqRange,TimeRange_DuringSent_Ind,TFCat_BaseNormed,MeanIndPowerTrialsCat);
        
        %%TimeRange_PosySent=[2.5 3.5];  %previuous selection: [2.5 4.5]; %Second
        MeanInducedPowerCat_PostSent(:,i)=FnFindInducedPowerinSelTime(FreqRange,TimeRange_PostSent_Ind,TFCat_BaseNormed,MeanIndPowerTrialsCat);
        
        FreqRange=[2 8]; %Hz
        %%TimeRange_PostMTB=[1 1.4];% Previuous selction: [1 1.5]; %Second
        MeanEvokedPowerCat_PostMTB(:,i)=FnFindEvokedPowerinSelTime(FreqRange,TimeRange_PostMTB_Evk,TFCat,EvokedPowerCat);
        
        %%TimeRange_PostSent=[2 2.4];%Second
        MeanEvokedPowerCat_PostSent(:,i)=FnFindEvokedPowerinSelTime(FreqRange,TimeRange_PostSent_Evk,TFCat,EvokedPowerCat);
        
        FreqRange=[8 12]; %Hz
        %%TimeRange_PreSent=[1.5 2]; %Second
        MeanEvokedPowerCat_PreSent(:,i)=FnFindEvokedPowerinSelTime(FreqRange,TimeRange_PreSent_Evk,TFCat,EvokedPowerCat);
        
    end
    if 1 %Use 0 for Grand Averaging
        %% Time-Frequency analysis
        if 0 %just test
            %*********************** Create S 14 (-5 dB SNR) **********************
            DataS14=FnCreateEpochedData(FileName,FinalEventName,FinalEventTimes,TimeRange,Fs,'S 14');
            %****************** Create S16 (+5 dB SNR)**************************
            DataS16=FnCreateEpochedData(FileName,FinalEventName,FinalEventTimes,TimeRange,Fs,'S 16');
        end
        %******************Time-Frequncy for  S14 ****************************
        %Update: Baseline normalization was added
        TFS14=FnTimeFreqAnalysis(DataS14,TFParam);
        TFS14.dimord='chan_freq_time';
        PowerTemp=FnInducedPower(TFS14.fourierspctrm);
        TFS14.powspctrm=squeeze(nanmean(PowerTemp,1));% I added power data for baseline norm.
        cfg=[];  cfg.baseline=BaselineWin; %previuous selection: [0.5 0.8]; %second
        cfg.parameter={'powspctrm'}; %fourierspctrm
        cfg.baselinetype=BaseLineNormtype; %or relchange %or dB
        %Condition-specefic besline
        %[TFS14_BaseNormed] = ft_freqbaseline(cfg, TFS14); %Baseline normalization
        %Averaged baseline
        % % data = 10*log10(data ./ meanVals);
        TFS14_BaseNormed=TFS14;
        if Norm
            TFS14_BaseNormed.powspctrm=10*log10( TFS14.powspctrm ./ meanVals);
        end
        %%InducedPowerCat=FnInducedPower(TFCat_BaseNormed.fourierspctrm);%Induced Power
        MeanIndPowerTrialsS14=  TFS14_BaseNormed.powspctrm;%Average over trilas

      



        
        if SelChanSwitch==1
            [SelLab,SelInd,~]=intersect(TFS14.label,SelChannLabels,'stable');%you should insert the desired channel labels, here
            MeanIndPowerTrialsChannS14_Sel(:,:,i)=squeeze(nanmean(MeanIndPowerTrialsS14(SelInd,:,:),1));%Average over channels
            [EvokedPowerS14,MeanEvokedTrialsChannS14_Sel(:,:,i)]=FnEvokedPower(TFS14.fourierspctrm(:,SelInd,:,:),Method);%Evoked power (for evoked power (plv), I did'nt apply baseline normalization, it does not make sense
        end
        MeanIndPowerTrialsChannS14(:,:,i)=squeeze(nanmean(MeanIndPowerTrialsS14,1));%Average over channels
        [EvokedPowerS14,MeanEvokedTrialsChannS14(:,:,i)]=FnEvokedPower(TFS14.fourierspctrm,Method);%Evoked power (for evoked power (plv), I did'nt apply baseline normalization, it does not make sense
        
        %         MeanIndPowerTrialsChannS14(:,:,i)=squeeze(nanmean(MeanIndPowerTrialsS14,1));%Average over channels
        %         [EvokedPowerS14,MeanEvokedTrialsChannS14(:,:,i)]=FnEvokedPower(TFS14.fourierspctrm,Method);%Evoked power (for evoked power (plv), I did'nt apply baseline normalization, it does not make sense
        
        
        
        
        % %         TFS14=FnTimeFreqAnalysis(DataS14,TFParam);
        % %         TFS14.dimord='chan_freq_time';
        % %         PowerTemp=FnInducedPower(TFS14.fourierspctrm);%induced power
        % %         TFS14.powspctrm=squeeze(nanmean(PowerTemp),1);% I added power data for baseline norm.
        % %         cfg=[];  cfg.baseline=BaselineWin; %previuous selection: [0.5 0.8]; %second
        % %         cfg.parameter={'powspctrm'}; %fourierspctrm
        % %         cfg.baselinetype='db'; %or relchange %or absolute
        % %         [TFS14_BaseNormed] = ft_freqbaseline(cfg, TFS14); %Baseline normalization
        % %         %%InducedPowerCat=FnInducedPower(TFCat_BaseNormed.fourierspctrm);%Induced Power
        % %         MeanIndPowerTrialsS14=TFS14_BaseNormed.powspctrm;%Induced Power
        % %         %         MeanIndPowerTrialsS14=squeeze(nanmean(InducedPowerS14,1));%Average over trilas
        % %         MeanIndPowerTrialsChannS14(:,:,i)=squeeze(nanmean(MeanIndPowerTrialsS14,1));%Average over channels
        % %         [EvokedPowerS14,MeanEvokedTrialsChannS14(:,:,i)]=FnEvokedPower(TFS14.fourierspctrm,Method);%Evoked power (for evoked power (plv), I did'nt apply baseline normalization, it does not make sense
        %******************Time-Frequncy for  S14 ****************************
        %Update: Baseline normalization was added
        
        TFS16=FnTimeFreqAnalysis(DataS16,TFParam);
        TFS16.dimord='chan_freq_time';
        PowerTemp=FnInducedPower(TFS16.fourierspctrm);
        TFS16.powspctrm=squeeze(nanmean(PowerTemp,1));% I added power data for baseline norm.
        cfg=[];  cfg.baseline=BaselineWin; %previuous selection: [0.5 0.8]; %second
        cfg.parameter={'powspctrm'}; %fourierspctrm
        cfg.baselinetype=BaseLineNormtype; %or relchange %or dB
        %Condition-specefic besline
        %[TFS16_BaseNormed] = ft_freqbaseline(cfg, TFS16); %Baseline normalization
        %Averaged baseline
        TFS16_BaseNormed=TFS16;
        if Norm
            TFS16_BaseNormed.powspctrm=10*log10( TFS16.powspctrm ./ meanVals);
        end
        %%InducedPowerCat=FnInducedPower(TFCat_BaseNormed.fourierspctrm);%Induced Power
        MeanIndPowerTrialsS16=  TFS16_BaseNormed.powspctrm;%Average over trilas
        %%%% Channel Selection
        if SelChanSwitch==1
            [SelLab,SelInd,~]=intersect(TFS16.label,SelChannLabels,'stable');%you should insert the desired channel labels, here
            MeanIndPowerTrialsChannS16_Sel(:,:,i)=squeeze(nanmean(MeanIndPowerTrialsS16(SelInd,:,:),1));%Average over channels
            [EvokedPowerS16,MeanEvokedTrialsChannS16_Sel(:,:,i)]=FnEvokedPower(TFS16.fourierspctrm(:,SelInd,:,:),Method);%Evoked power (for evoked power (plv), I did'nt apply baseline normalization, it does not make sense
        end
            MeanIndPowerTrialsChannS16(:,:,i)=squeeze(nanmean(MeanIndPowerTrialsS16,1));%Average over channels
            [EvokedPowerS16,MeanEvokedTrialsChannS16(:,:,i)]=FnEvokedPower(TFS16.fourierspctrm,Method);%Evoked power (for evoked power (plv), I did'nt apply baseline normalization, it does not make sense
        
        
        
        % %         MeanIndPowerTrialsChannS16(:,:,i)=squeeze(nanmean(MeanIndPowerTrialsS16,1));%Average over channels
        % %         [EvokedPowerS16,MeanEvokedTrialsChannS16(:,:,i)]=FnEvokedPower(TFS16.fourierspctrm,Method);%Evoked power (for evoked power (plv), I did'nt apply baseline normalization, it does not make sense
        
        
        % %         TFS16=FnTimeFreqAnalysis(DataS16,TFParam);
        % %         TFS16.dimord='chan_freq_time';
        % %         PowerTemp=FnInducedPower(TFS16.fourierspctrm);%induced power
        % %         TFS16.powspctrm=squeeze(nanmean(PowerTemp),1);% I added power data for baseline norm.
        % %         cfg=[];  cfg.baseline=BaselineWin; %previuous selection: [0.5 0.8]; %second
        % %         cfg.parameter={'powspctrm'}; %fourierspctrm
        % %         cfg.baselinetype='db'; %or relchange %or absolute
        % %         [TFS16_BaseNormed] = ft_freqbaseline(cfg, TFS16); %Baseline normalization
        % %         %%InducedPowerCat=FnInducedPower(TFCat_BaseNormed.fourierspctrm);%Induced Power
        % %         MeanIndPowerTrialsS16=TFS16_BaseNormed.powspctrm;%Induced Power
        % %         %         MeanIndPowerTrialsS16=squeeze(nanmean(InducedPowerS16,1));%Average over trilas
        % %         MeanIndPowerTrialsChannS16(:,:,i)=squeeze(nanmean(MeanIndPowerTrialsS16,1));%Average over channels
        % %         [EvokedPowerS16,MeanEvokedTrialsChannS16(:,:,i)]=FnEvokedPower(TFS16.fourierspctrm,Method);%Evoked power (for evoked power (plv), I did'nt apply baseline normalization, it does not make sense
        %%
        if 0%Use 1 if you want to run the code based on the preliminary results (26/11/2020)
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
            %         if 0
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
            caxis([min(MeanPowerS14(:)) max(MeanPowerS14(:))]) %-max(MeanPowerS14(:))
            shading interp;
            xlabel('Time (sec)');ylabel('Frequency (Hz)');zlabel('Power spectrum');
            view(0,-90)
            h=colorbar;
            colormap('jet');
            ylabel(h, 'Power averaged across all channels for the whole trial duration')
            title(['+5 dB SNR, S16-Subject-' num2str(i)])
            %         end
            %% Evoked Power
            %**************************** Evoked Power S14 ************************
            [EvokedPowerS14,MeanEvokedPowerS14]=FnEvokedPower(TFS14.fourierspctrm,Method);
            %**************************** Evoked Power S16 ************************
            [EvokedPowerS16,MeanEvokedPowerS16]=FnEvokedPower(TFS16.fourierspctrm,Method);
            %*******************************Plot S14 evoked power******************
            %         if 0
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
        end
        % ***************************** pre-stimulus induced alpha-S 14*************
        %Upadate: useing baseline normalized TF
        FreqRange=[8 12]; %Hz (alpha range)
        %%TimeRange_PreSent= [1.5 2]; %previuous selection: [1.4 1.9]; %Second
        MeanInducedPowerS14_PreSent(:,i)=FnFindInducedPowerinSelTime(FreqRange,TimeRange_PreSent_Ind,TFS14_BaseNormed,MeanIndPowerTrialsS14);
        %Update: coherent power (15/12/2024)
        %Add coherent power (Eig on Coh)
       % Perform eigendecomposition

       % [XYw,Cmat,Ctot,Cvec,Cent,Sdiag]=FnEigCrossSpectrum(Xw); % Calculate cross-spectra
       TF142Cut=[];TF142Cut=TFS14_BaseNormed.fourierspctrm;
       XwS14=[];XwS14=FnFindInducedPowerinSelTime(FreqRange,TimeRange_PreSent_Ind,TFS14_BaseNormed,TF142Cut);
       XwS14=permute(XwS14,[2 1 3]);
       [XYw_S14,Cmat_S14,Ctot_S14,Cvec_S14,Cent_S14,Sdiag_S14]=FnEigCrossSpectrum(XwS14); % Calculate cross-spectra
        PreSentCoherentAlphaPower_S14(:,i)=nanmean(abs(Cvec_S14),1);
        PreSentSpectralAlphaPower_S14(:,i)=nanmean(Sdiag_S14,1);




   
        
        %%TimeRange_DuringSent=[2.5 3.5];  %previuous selection: [2.5 4.5]; %Second
        MeanInducedPowerS14_DuringSent(:,i)=FnFindInducedPowerinSelTime(FreqRange,TimeRange_DuringSent_Ind,TFS14_BaseNormed,MeanIndPowerTrialsS14);
        %Update: coherent power (15/12/2024)
        %Add coherent power (Eig on Coh)
        % Perform eigendecomposition

        % [XYw,Cmat,Ctot,Cvec,Cent,Sdiag]=FnEigCrossSpectrum(Xw); % Calculate cross-spectra
        TF142Cut=[];TF142Cut=TFS14_BaseNormed.fourierspctrm;
        XwS14=[];XwS14=FnFindInducedPowerinSelTime(FreqRange,TimeRange_DuringSent_Ind,TFS14_BaseNormed,TF142Cut);
        XwS14=permute(XwS14,[2 1 3]);
        [XYw_S14,Cmat_S14,Ctot_S14,Cvec_S14,Cent_S14,Sdiag_S14]=FnEigCrossSpectrum(XwS14); % Calculate cross-spectra
        DurSentCoherentAlphaPower_S14(:,i)=nanmean(abs(Cvec_S14),1);
        DurSentSpectralAlphaPower_S14(:,i)=nanmean(Sdiag_S14,1);


        %%TimeRange_PostSent=[4.5 6]; %Second
        MeanInducedPowerS14_PostSent(:,i)=FnFindInducedPowerinSelTime(FreqRange,TimeRange_PostSent_Ind,TFS14_BaseNormed,MeanIndPowerTrialsS14);%This power probably is not necessery to check
        %Update: coherent power (15/12/2024)
        %Add coherent power (Eig on Coh)
        % Perform eigendecomposition

        % [XYw,Cmat,Ctot,Cvec,Cent,Sdiag]=FnEigCrossSpectrum(Xw); % Calculate cross-spectra
        TF142Cut=[];TF142Cut=TFS14_BaseNormed.fourierspctrm;
        XwS14=[];XwS14=FnFindInducedPowerinSelTime(FreqRange,TimeRange_PostSent_Ind,TFS14_BaseNormed,TF142Cut);
        XwS14=permute(XwS14,[2 1 3]);
        [XYw_S14,Cmat_S14,Ctot_S14,Cvec_S14,Cent_S14,Sdiag_S14]=FnEigCrossSpectrum(XwS14); % Calculate cross-spectra
        PostSentCoherentAlphaPower_S14(:,i)=nanmean(abs(Cvec_S14),1);
        PostSentSpectralAlphaPower_S14(:,i)=nanmean(Sdiag_S14,1);
        % ***************************** Plot pre-stimulus induced alpha- S 14************
        if 0 %Use 1 for visulaization
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
        end
        % ***************************** pre-stimulus induced alpha-S 16*************
        %Upadate: useing baseline normalized TF
        FreqRange=[8 12]; %Hz (alpha range)
        %%TimeRange_PreSent= [1.5 2.0]; %previuous selection: [1.4 1.9]; %Second
        MeanInducedPowerS16_PreSent(:,i)=FnFindInducedPowerinSelTime(FreqRange,TimeRange_PreSent_Ind,TFS16_BaseNormed,MeanIndPowerTrialsS16);
        
        %%TimeRange_DuringSent=[2.5 3.5];  %previuous selection: [2.5 4.5]; %Second
        MeanInducedPowerS16_DuringSent(:,i)=FnFindInducedPowerinSelTime(FreqRange,TimeRange_DuringSent_Ind,TFS16_BaseNormed,MeanIndPowerTrialsS16);
        
        
        %         TimeRange_DuringSent=[4.5 6]; %Second
        MeanInducedPowerS16_PostSent(:,i)=FnFindInducedPowerinSelTime(FreqRange,TimeRange_PostSent_Ind,TFS16_BaseNormed,MeanIndPowerTrialsS16);%This power probably is not necessery to check
        % ***************************** Plot pre-stimulus induced alpha- S 14************
        if 0 %Use 1 for visulaization
            FnTopoPlotPower(MeanInducedPowerS16_PreSent,EEGChanLoc,0,max(MeanInducedPowerS14_PreSent),Type);
            title(['+5 dB SNR, S16-Subject-' num2str(i)  '-[' num2str(TimeRange_PreSent(1)) '-'...
                num2str(TimeRange_PreSent(2)) '] sec-[' num2str(FreqRange(1)) '-' num2str(FreqRange(2)) '] Hz']);
            
            FnTopoPlotPower(MeanInducedPowerS16_DuringSent,EEGChanLoc,0,max(MeanInducedPowerS14_DuringSent),Type);
            title(['+5 dB SNR, S16-Subject-' num2str(i)  '-[' num2str(TimeRange_DuringSent(1)) '-'...
                num2str(TimeRange_DuringSent(2)) '] sec-[' num2str(FreqRange(1)) '-' num2str(FreqRange(2)) '] Hz']);
            
            FnTopoPlotPower(MeanInducedPowerS16_PostSent,EEGChanLoc,0,max(MeanInducedPowerS14_PostSent),Type);
            title(['+5 dB SNR, S16-Subject-' num2str(i)  '-[' num2str(TimeRange_PostSent(1)) '-'...
                num2str(TimeRange_PostSent(2)) '] sec-[' num2str(FreqRange(1)) '-' num2str(FreqRange(2)) '] Hz']);
        end
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
        %%TimeRange_PostMTB=[1 1.4];% Previuous selction: [1 1.5]; %Second
        MeanEvokedPowerS14_PostMTB(:,i)=FnFindEvokedPowerinSelTime(FreqRange,TimeRange_PostMTB_Evk,TFS14,EvokedPowerS14);
        
        %%TimeRange_PostSent=[2 2.4]; %Second
        MeanEvokedPowerS14_PostSent(:,i)=FnFindEvokedPowerinSelTime(FreqRange,TimeRange_PostSent_Evk,TFS14,EvokedPowerS14);
        % ***************************** Post-MTB evoked  LF power S16*************
        FreqRange=[2 8]; %Hz
        
        FreqRange_alpha=[8 12]; %Hz
        %%TimeRange_PreMTB=[1.5 2];%Second
        MeanEvokedPowerS14_PreSent(:,i)=FnFindEvokedPowerinSelTime(FreqRange_alpha,TimeRange_PreSent_Evk,TFS14,EvokedPowerS14);
    
        
        %% TimeRange_PostMTB=[1 1.4];% Previuous selction: [1 1.5]; %Second
        MeanEvokedPowerS16_PostMTB(:,i)=FnFindEvokedPowerinSelTime(FreqRange,TimeRange_PostMTB_Evk,TFS16_BaseNormed,EvokedPowerS16);
        
        %%TimeRange_PostSent=[2 2.4]; %Second
        MeanEvokedPowerS16_PostSent(:,i)=FnFindEvokedPowerinSelTime(FreqRange,TimeRange_PostSent_Evk,TFS16_BaseNormed,EvokedPowerS16); %This power probably is not necessery to check
 
        %%TimeRange_PreMTB=[1.5 2];%Second
        MeanEvokedPowerS16_PreSent(:,i)=FnFindEvokedPowerinSelTime(FreqRange_alpha,TimeRange_PreSent_Evk,TFS16,EvokedPowerS16);
        % ***************************** Plot evoked alpha- S 14************
        if 0 %use 1 for visulaization and saving the results based on prelimanary results
            
        
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
            %         end
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
                MeanEvokedPowerS16_PostMTB MeanEvokedPowerS16_PostSent...
                MeanEvokedPowerS14_PreSent...
                MeanEvokedPowerS16_PreSent];
            
            VarNames={'Channles','IndS14PreSent','IndS14DurSent','IndS14PostSent',...
                'IndS16PreSent','IndS16DurSent','IndS16PostSent'...
                'EvS14PostMTB','EvS14PostSent',...
                'EvS16PostMTB','EvS16PostSent', ...
                'Ev14PreSent', 'Ev16PreSent'};%Variable names
            Table = table(TFS14.label,MeanInducedPowerS14_PreSent,MeanInducedPowerS14_DuringSent,MeanInducedPowerS14_PostSent, ...
                MeanInducedPowerS16_PreSent,MeanInducedPowerS16_DuringSent,MeanInducedPowerS16_PostSent,...
                MeanEvokedPowerS14_PostMTB,MeanEvokedPowerS14_PostSent,...
                MeanEvokedPowerS16_PostMTB,MeanEvokedPowerS16_PostSent, MeanEvokedPowerS14_PreSent, MeanEvokedPowerS16_PreSent);%put in a table
            Table.Properties.VariableNames=VarNames;%assign the variable names
            writetable(Table,[EDFDir(i).name(1:end-4) '_TimeFreqPower' '.txt'],'Delimiter',' '); %write
        end
        
        %% Writing the results (Update: 12/03/2021)
        Table=[];
        %VarNames={'SubjectID' ,'Group' ,'Condition' ,'Channel','PreSentIndPow','DurSentIndPow','PostSentIndPow','PostMTBEvkPow','PostSentEvkPow'};
        CurrentName=EDFDir(i).name;
        ULIndCurrent=strfind(CurrentName,'_');
        GroupandID=CurrentName(ULIndCurrent(1)+1:ULIndCurrent(2)-1);
        ChanNum=length(MeanInducedPowerS14_PreSent(:,i));
        Table.SubjectID=cellstr(repmat(GroupandID(3:4),[ChanNum*2,1]));
        Table.Group=cellstr(repmat(GroupandID(1:2),[ChanNum*2,1]));
        ChanLab=OurLayout.label';
        Table.Channel=repmat(ChanLab,[2,1]);
        Table.Condition=cell(ChanNum*2,1);
        Table.Condition(1:ChanNum)=cellstr(repmat('S14',[ChanNum,1]));
        Table.Condition(ChanNum+1:end)=cellstr(repmat('S16',[ChanNum,1]));
        Table.PreSentIndPow=repmat(MeanInducedPowerS14_PreSent(:,i),[2,1]);
        Table.PreSentIndPow(ChanNum+1:end)=MeanInducedPowerS16_PreSent(:,i);
        
        Table.DurSentIndPow=repmat(MeanInducedPowerS14_DuringSent(:,i),[2,1]);
        Table.DurSentIndPow(ChanNum+1:end)=MeanInducedPowerS16_DuringSent(:,i);
        
        Table.PostSentIndPow=repmat(MeanInducedPowerS14_PostSent(:,i),[2,1]);
        Table.PostSentIndPow(ChanNum+1:end)=MeanInducedPowerS16_PostSent(:,i);
        
        Table.PostMTBEvkPow=repmat(MeanEvokedPowerS14_PostMTB(:,i),[2,1]);
        Table.PostMTBEvkPow(ChanNum+1:end)=MeanEvokedPowerS16_PostMTB(:,i);
        
        Table.PostSentEvkPow=repmat(MeanEvokedPowerS14_PostSent(:,i),[2,1]);
        Table.PostSentEvkPow(ChanNum+1:end)=MeanEvokedPowerS16_PostSent(:,i);
        %% Added these codes after paper revision (pre-sent...)
        Table.PreSentEvkPow=repmat(MeanEvokedPowerS14_PreSent(:,i),[2,1]);
        Table.PreSentEvkPow(ChanNum+1:end)=MeanEvokedPowerS16_PreSent(:,i);
    end
    CatTable=[CatTable;struct2table(Table)];
    
end
writetable(CatTable,['TimeFreqPower_Features.txt'],'Delimiter',' '); %write
writetable(CatTable,['TimeFreqPower_Features.xls']); %write

if 1
    %%  Grand-average time-frequaency represneation across trial, conditions, subjects, and electrodes
    GrandAvgInducedTF=nanmean(MeanIndPowerTrialsChannCat,3);
    GrandAvgEvokedTF=nanmean(MeanEvokedTrialsChannCat,3);
    
    figure,surf(TFParam.toi,TFParam.foi,GrandAvgInducedTF);
    caxis([min(min(GrandAvgInducedTF(:)),min(GrandAvgInducedTF(:))) max(max(GrandAvgInducedTF(:)),max(GrandAvgInducedTF(:)))])
    shading interp;
    xlabel('Time (sec)');ylabel('Frequency (Hz)');zlabel('Power spectrum');
    view(0,-90)
    h=colorbar;
    colormap('jet');
    ylabel(h, ['Induced Grand-average TF across trial, subjects, conditions and electrodes']);
    title('Merged Conditions-Induced')
    
    figure,surf(TFParam.toi,TFParam.foi,GrandAvgEvokedTF);
    caxis([min(min(GrandAvgEvokedTF(:)),min(GrandAvgEvokedTF(:))) max(max(GrandAvgEvokedTF(:)),max(GrandAvgEvokedTF(:)))])
    shading interp;
    xlabel('Time (sec)');ylabel('Frequency (Hz)');zlabel('Power spectrum');
    view(0,-90)
    h=colorbar;
    colormap('jet');
    ylabel(h, ['Evoked Grand-average TF across trial, subjects, conditions and electrodes']);
    title('Merged Conditions-Evoked')
    
    %%%% Conditions
    %S14
    GrandAvgInducedTFS14=nanmean(MeanIndPowerTrialsChannS14,3);
    GrandAvgEvokedTFS14=nanmean(MeanEvokedTrialsChannS14,3);
    %S16
    GrandAvgInducedTFS16=nanmean(MeanIndPowerTrialsChannS16,3);
    GrandAvgEvokedTFS16=nanmean(MeanEvokedTrialsChannS16,3);
    
    figure,surf(TFParam.toi,TFParam.foi,GrandAvgInducedTFS14);
    caxis([min(min(GrandAvgInducedTF(:)),min(GrandAvgInducedTF(:))) max(max(GrandAvgInducedTF(:)),max(GrandAvgInducedTF(:)))])
    shading interp;
    xlabel('Time (sec)');ylabel('Frequency (Hz)');zlabel('Power spectrum');
    view(0,-90)
    h=colorbar;
    colormap('jet');
    ylabel(h, ['Induced S14 TF across trial, subjects and electrodes']);
    title('-5 dB-Induced Power');
    
    
    figure,surf(TFParam.toi,TFParam.foi,GrandAvgEvokedTFS14);
    caxis([min(min(GrandAvgEvokedTF(:)),min(GrandAvgEvokedTF(:))) max(max(GrandAvgEvokedTF(:)),max(GrandAvgEvokedTF(:)))])
    shading interp;
    xlabel('Time (sec)');ylabel('Frequency (Hz)');zlabel('Power spectrum');
    view(0,-90)
    h=colorbar;
    colormap('jet');
    ylabel(h, ['Evoked S16 TF across trial, subjects and electrodes']);
    title('-5 dB-Evoked Power');
    
    
    figure,surf(TFParam.toi,TFParam.foi,GrandAvgInducedTFS16);
    caxis([min(min(GrandAvgInducedTF(:)),min(GrandAvgInducedTF(:))) max(max(GrandAvgInducedTF(:)),max(GrandAvgInducedTF(:)))])
    shading interp;
    xlabel('Time (sec)');ylabel('Frequency (Hz)');zlabel('Power spectrum');
    view(0,-90)
    h=colorbar;
    colormap('jet');
    ylabel(h, ['Induced S16 TF across trial, subjects and electrodes']);
    title('+5 dB-Induced Power');
    
    
    figure,surf(TFParam.toi,TFParam.foi,GrandAvgEvokedTFS16);
    caxis([min(min(GrandAvgEvokedTF(:)),min(GrandAvgEvokedTF(:))) max(max(GrandAvgEvokedTF(:)),max(GrandAvgEvokedTF(:)))])
    shading interp;
    xlabel('Time (sec)');ylabel('Frequency (Hz)');zlabel('Power spectrum');
    view(0,-90)
    h=colorbar;
    colormap('jet');
    ylabel(h, ['Evoked S16 TF across trial, subjects and electrodes']);
    title('+5 dB-Evoked Power');
    %% SelChan plot
    if SelChanSwitch==1
        GrandAvgInducedTFS14_Sel=nanmean(MeanIndPowerTrialsChannS14_Sel,3);
        GrandAvgEvokedTFS14_Sel=nanmean(MeanEvokedTrialsChannS14_Sel,3);
        %S16
        GrandAvgInducedTFS16_Sel=nanmean(MeanIndPowerTrialsChannS16_Sel,3);
        GrandAvgEvokedTFS16_Sel=nanmean(MeanEvokedTrialsChannS16_Sel,3);
        
        figure,surf(TFParam.toi,TFParam.foi,GrandAvgInducedTFS14_Sel);
        caxis([min(min(GrandAvgInducedTF(:)),min(GrandAvgInducedTF(:))) max(max(GrandAvgInducedTF(:)),max(GrandAvgInducedTF(:)))])
        shading interp;
        xlabel('Time (sec)');ylabel('Frequency (Hz)');zlabel('Power spectrum');
        view(0,-90)
        h=colorbar;
        colormap('jet');
        ylabel(h, ['Induced S14 TF across trial, subjects and electrodes']);
        title('-5 dB-Induced PowerSelcted channels');
        
        
        figure,surf(TFParam.toi,TFParam.foi,GrandAvgEvokedTFS14_Sel);
        caxis([min(min(GrandAvgEvokedTF(:)),min(GrandAvgEvokedTF(:))) max(max(GrandAvgEvokedTF(:)),max(GrandAvgEvokedTF(:)))])
        shading interp;
        xlabel('Time (sec)');ylabel('Frequency (Hz)');zlabel('Power spectrum');
        view(0,-90)
        h=colorbar;
        colormap('jet');
        ylabel(h, ['Evoked S16 TF across trial, subjects and electrodes']);
        title('-5 dB-Evoked PowerSelcted channels');
        
        
        figure,surf(TFParam.toi,TFParam.foi,GrandAvgInducedTFS16_Sel);
        caxis([min(min(GrandAvgInducedTF(:)),min(GrandAvgInducedTF(:))) max(max(GrandAvgInducedTF(:)),max(GrandAvgInducedTF(:)))])
        shading interp;
        xlabel('Time (sec)');ylabel('Frequency (Hz)');zlabel('Power spectrum');
        view(0,-90)
        h=colorbar;
        colormap('jet');
        ylabel(h, ['Induced S16 TF across trial, subjects and electrodes']);
        title('+5 dB-Induced Power- Selcted channels');
        
        
        figure,surf(TFParam.toi,TFParam.foi,GrandAvgEvokedTFS16_Sel);
        caxis([min(min(GrandAvgEvokedTF(:)),min(GrandAvgEvokedTF(:))) max(max(GrandAvgEvokedTF(:)),max(GrandAvgEvokedTF(:)))])
        shading interp;
        xlabel('Time (sec)');ylabel('Frequency (Hz)');zlabel('Power spectrum');
        view(0,-90)
        h=colorbar;
        colormap('jet');
        ylabel(h, ['Evoked S16 TF across trial, subjects and electrodes']);
        title('+5 dB-Evoked PowerSelcted channels');
    end
    %%  Grand-average topographic maps across trial, conditions, subjects, desired frequency and time
    GrandAvg_InducedPowerCat_PreSent=nanmean(MeanInducedPowerCat_PreSent,2);
    GrandAvg_InducedPowerCat_DuringSent=nanmean(MeanInducedPowerCat_DuringSent,2);
    
    GrandAvg_InducedPowerCat_PostSent=nanmean(MeanInducedPowerCat_PostSent,2);
    
    GrandAvg_EvokedPowerCat_PostMTB=nanmean(MeanEvokedPowerCat_PostMTB,2);
    GrandAvg_EvokedPowerCat_PostSent=nanmean(MeanEvokedPowerCat_PostSent,2);
    
    %topoplot for gg-av
    Type='GrandAvg Induced Power PreSent';
    FnTopoPlotPower(GrandAvg_InducedPowerCat_PreSent,EEGChanLoc,-3,3,Type);
    title(['Induced Power PreSent'   '-[' num2str(TimeRange_PreSent_Ind(1)) '-'...
        num2str(TimeRange_PreSent_Ind(2)) '] sec-[' num2str(8) '-' num2str(12) '] Hz']);
    
    
    Type='GrandAvg Induced Power DuringSent';
    FnTopoPlotPower(GrandAvg_InducedPowerCat_DuringSent,EEGChanLoc,-3,3,Type);
    title(['Induced Power DuringSent'  '-[' num2str(TimeRange_DuringSent_Ind(1)) '-'...
        num2str(TimeRange_DuringSent_Ind(2)) '] sec-[' num2str(8) '-' num2str(12) '] Hz']);
    
    
    
    Type='GrandAvg Induced Power PostSent';
    FnTopoPlotPower(GrandAvg_InducedPowerCat_PostSent,EEGChanLoc,0,1,Type);
    title(['Induced Power DuringSent'  '-[' num2str(TimeRange_PostSent_Ind(1)) '-'...
        num2str(TimeRange_PostSent_Ind(2)) '] sec-[' num2str(8) '-' num2str(12) '] Hz']);
    %%Evoked
    Type='GrandAvg Evoked Power PostMTB';
    FnTopoPlotPower(GrandAvg_EvokedPowerCat_PostMTB,EEGChanLoc,0,.5,Type);
    title(['Evoked Power PostMTB'   '-[' num2str(TimeRange_PostMTB_Evk(1)) '-'...
        num2str(TimeRange_PostMTB_Evk(2)) '] sec-[' num2str(2) '-' num2str(8) '] Hz']);
    
    Type='GrandAvg Evoked Power PostSentence';
    FnTopoPlotPower(GrandAvg_EvokedPowerCat_PostSent,EEGChanLoc,0,.5,Type);
    title(['Evoked Power PostSentence'   '-[' num2str(TimeRange_PostSent_Evk(1)) '-'...
        num2str(TimeRange_PostSent_Evk(2)) '] sec-[' num2str(2) '-' num2str(8) '] Hz']);
    %%%Conditions
    %% S14
    %%%%% Induced Power
    GrandAvg_InducedPowerS14_PreSent=nanmean(MeanInducedPowerS14_PreSent,2);
    GrandAvg_InducedPowerS14_DuringSent=nanmean(MeanInducedPowerS14_DuringSent,2);
    GrandAvg_InducedPowerS14_PostSent=nanmean(MeanInducedPowerS14_PostSent,2);

    GrandAvg_AlphaCohPowerS14_PreSent=nanmean(PreSentCoherentAlphaPower_S14,2);
    GrandAvg_AlphaCohPowerS14_DuringSent=nanmean(DurSentCoherentAlphaPower_S14,2);
    GrandAvg_AlphaCohPowerS14_PostSent=nanmean(PostSentCoherentAlphaPower_S14,2);

    GrandAvg_AlphaSpecPowerS14_PreSent=nanmean(PreSentSpectralAlphaPower_S14,2);
    GrandAvg_AlphaSpecPowerS14_DuringSent=nanmean(DurSentSpectralAlphaPower_S14,2);
    GrandAvg_AlphaSpecPowerS14_PostSent=nanmean(DurSentSpectralAlphaPower_S14,2);

    %%%%%Evoked Power
    GrandAvg_EvokedPowerS14_PostMTB=nanmean(MeanEvokedPowerS14_PostMTB,2);
    GrandAvg_EvokedPowerS14_PostSent=nanmean(MeanEvokedPowerS14_PostSent,2);
    
    
    %topoplot
    Type='GrandAvg Induced Power PreSent';
    FnTopoPlotPower(GrandAvg_InducedPowerS14_PreSent,EEGChanLoc,min(GrandAvg_InducedPowerS14_PreSent),max(GrandAvg_InducedPowerS14_PreSent),Type);%0,2
    title(['-5dB-Induced Power PreSent'   '-[' num2str(TimeRange_PreSent_Ind(1)) '-'...
        num2str(TimeRange_PreSent_Ind(2)) '] sec-[' num2str(8) '-' num2str(12) '] Hz']);
    
    
    Type='GrandAvg Induced Power DuringSent';
    FnTopoPlotPower(GrandAvg_InducedPowerS14_DuringSent,EEGChanLoc,min(GrandAvg_InducedPowerS14_DuringSent),max(GrandAvg_InducedPowerS14_DuringSent),Type);%-3,3
    title(['-5dB-Induced Power DuringSent'  '-[' num2str(TimeRange_DuringSent_Ind(1)) '-'...
        num2str(TimeRange_DuringSent_Ind(2)) '] sec-[' num2str(8) '-' num2str(12) '] Hz']);
    
    
    
    Type='GrandAvg Induced Power PostSent';
    FnTopoPlotPower(GrandAvg_InducedPowerS14_PostSent,EEGChanLoc,min(GrandAvg_InducedPowerS14_DuringSent),max(GrandAvg_InducedPowerS14_DuringSent),Type);%0,1
    title(['-5dB-Induced Power DuringSent'  '-[' num2str(TimeRange_PostSent_Ind(1)) '-'...
        num2str(TimeRange_PostSent_Ind(2)) '] sec-[' num2str(8) '-' num2str(12) '] Hz']);
    %Coh

    Type='GrandAvg Alpha-Coh Power PreSent';
    FnTopoPlotPower(GrandAvg_AlphaCohPowerS14_PreSent,EEGChanLoc,min(GrandAvg_AlphaCohPowerS14_PreSent),max(GrandAvg_AlphaCohPowerS14_PreSent),Type);
    title(['-5dB-Induced Power PreSent'   '-[' num2str(TimeRange_PreSent_Ind(1)) '-'...
        num2str(TimeRange_PreSent_Ind(2)) '] sec-[' num2str(8) '-' num2str(12) '] Hz']);


    Type='GrandAvg Alpha-Coh Power DuringSent';
    FnTopoPlotPower(GrandAvg_AlphaCohPowerS14_DuringSent,EEGChanLoc,min(GrandAvg_AlphaCohPowerS14_PreSent),max(GrandAvg_AlphaCohPowerS14_PreSent),Type);
    title(['-5dB-Induced Power DuringSent'  '-[' num2str(TimeRange_DuringSent_Ind(1)) '-'...
        num2str(TimeRange_DuringSent_Ind(2)) '] sec-[' num2str(8) '-' num2str(12) '] Hz']);



    Type='GrandAvg Alpha-Coh Power PostSent';
    FnTopoPlotPower(GrandAvg_AlphaCohPowerS14_PostSent,EEGChanLoc,min(GrandAvg_AlphaCohPowerS14_PreSent),max(GrandAvg_AlphaCohPowerS14_PreSent),Type);
    title(['-5dB-Induced Power PostSent'  '-[' num2str(TimeRange_PostSent_Ind(1)) '-'...
        num2str(TimeRange_PostSent_Ind(2)) '] sec-[' num2str(8) '-' num2str(12) '] Hz']);


    %SDiag

    Type='GrandAvg Alpha-Spec Power PreSent';
    FnTopoPlotPower(GrandAvg_AlphaSpecPowerS14_PreSent,EEGChanLoc,min(GrandAvg_AlphaSpecPowerS14_PreSent),max(GrandAvg_AlphaSpecPowerS14_PreSent),Type);
    title(['-5dB-Induced Power PreSent'   '-[' num2str(TimeRange_PreSent_Ind(1)) '-'...
        num2str(TimeRange_PreSent_Ind(2)) '] sec-[' num2str(8) '-' num2str(12) '] Hz']);


    Type='GrandAvg Alpha-Spec Power DuringSent';
    FnTopoPlotPower(GrandAvg_AlphaSpecPowerS14_DuringSent,EEGChanLoc,min(GrandAvg_AlphaSpecPowerS14_DuringSent),max(GrandAvg_AlphaSpecPowerS14_DuringSent),Type);
    title(['-5dB-Induced Power DuringSent'  '-[' num2str(TimeRange_DuringSent_Ind(1)) '-'...
        num2str(TimeRange_DuringSent_Ind(2)) '] sec-[' num2str(8) '-' num2str(12) '] Hz']);



    Type='GrandAvg Alpha-Spec Power PostSent';
    FnTopoPlotPower(GrandAvg_AlphaSpecPowerS14_PostSent,EEGChanLoc,min(GrandAvg_AlphaSpecPowerS14_PostSent),max(GrandAvg_AlphaSpecPowerS14_PostSent),Type);
    title(['-5dB-Induced Power PostSent'  '-[' num2str(TimeRange_PostSent_Ind(1)) '-'...
        num2str(TimeRange_PostSent_Ind(2)) '] sec-[' num2str(8) '-' num2str(12) '] Hz']);


    %%Evoked
    Type='GrandAvg Evoked Power PostMTB';
    FnTopoPlotPower(GrandAvg_EvokedPowerS14_PostMTB,EEGChanLoc,0,.5,Type);
    title(['-5dB-Evoked Power PostMTB'   '-[' num2str(TimeRange_PostMTB_Evk(1)) '-'...
        num2str(TimeRange_PostMTB_Evk(2)) '] sec-[' num2str(2) '-' num2str(8) '] Hz']);
    
    Type='GrandAvg Evoked Power PostSentence';
    FnTopoPlotPower(GrandAvg_EvokedPowerS14_PostSent,EEGChanLoc,0,.5,Type);
    title(['-5dB-Evoked Power PostSentence'   '-[' num2str(TimeRange_PostSent_Evk(1)) '-'...
        num2str(TimeRange_PostSent_Evk(2)) '] sec-[' num2str(2) '-' num2str(8) '] Hz']);
    
    %% S16
    %%%%% Induced Power
    GrandAvg_InducedPowerS16_PreSent=nanmean(MeanInducedPowerS16_PreSent,2);
    GrandAvg_InducedPowerS16_DuringSent=nanmean(MeanInducedPowerS16_DuringSent,2);
    GrandAvg_InducedPowerS16_PostSent=nanmean(MeanInducedPowerS16_PostSent,2);
    %%%%% Evoked Power
    GrandAvg_EvokedPowerS16_PostMTB=nanmean(MeanEvokedPowerS16_PostMTB,2);
    GrandAvg_EvokedPowerS16_PostSent=nanmean(MeanEvokedPowerS16_PostSent,2);
    %topoplot
    Type='GrandAvg Induced Power PreSent';
    FnTopoPlotPower(GrandAvg_InducedPowerS16_PreSent,EEGChanLoc,0,2,Type);
    title(['+5dB-Induced Power PreSent'   '-[' num2str(TimeRange_PreSent_Ind(1)) '-'...
        num2str(TimeRange_PreSent_Ind(2)) '] sec-[' num2str(8) '-' num2str(12) '] Hz']);
    
    
    Type='GrandAvg Induced Power DuringSent';
    FnTopoPlotPower(GrandAvg_InducedPowerS16_DuringSent,EEGChanLoc,0,2,Type);
    title(['+5dB-Induced Power DuringSent'  '-[' num2str(TimeRange_DuringSent_Ind(1)) '-'...
        num2str(TimeRange_DuringSent_Ind(2)) '] sec-[' num2str(8) '-' num2str(12) '] Hz']);
    
    
    
    Type='GrandAvg Induced Power PostSent';
    FnTopoPlotPower(GrandAvg_InducedPowerS16_PostSent,EEGChanLoc,0,2,Type);
    title(['+5dB-Induced Power DuringSent'  '-[' num2str(TimeRange_PostSent_Ind(1)) '-'...
        num2str(TimeRange_PostSent_Ind(2)) '] sec-[' num2str(8) '-' num2str(12) '] Hz']);
    %%Evoked
    Type='GrandAvg Evoked Power PostMTB';
    FnTopoPlotPower(GrandAvg_EvokedPowerS16_PostMTB,EEGChanLoc,0,.5,Type);
    title(['+5dB-Evoked Power PostMTB'   '-[' num2str(TimeRange_PostMTB_Evk(1)) '-'...
        num2str(TimeRange_PostMTB_Evk(2)) '] sec-[' num2str(2) '-' num2str(8) '] Hz']);
    
    Type='GrandAvg Evoked Power PostSentence';
    FnTopoPlotPower(GrandAvg_EvokedPowerS16_PostSent,EEGChanLoc,0,.5,Type);
    title(['+5dB-Evoked Power PostSentence'   '-[' num2str(TimeRange_PostSent_Evk(1)) '-'...
        num2str(TimeRange_PostSent_Evk(2)) '] sec-[' num2str(2) '-' num2str(8) '] Hz']);
    
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
end
%% output to txt --> excell
%% time windows : evoked [1 1.5] and [2 2.5]; indiced alpha [1.5 2] [2.5 4.5] [4.5 6]

