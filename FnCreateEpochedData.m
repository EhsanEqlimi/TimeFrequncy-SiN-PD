function Data=FnCreateEpochedData (FileName,FinalEventName,FinalEventTimes,TimeRange,Fs,EventName)
cfg=[];
% cfg.dataset=[FileName(1:end-4) '_ChannelAdded.edf'];
cfg.dataset=FileName;
cfg.continuous = 'yes';
cfg.channel= 'all';
S=FinalEventTimes(ismember(FinalEventName,EventName));
NTrial=length(S);
cfg.trl =[S+(TimeRange(1)*Fs),S+(TimeRange(2)*Fs), zeros(NTrial,1)];
Data=ft_preprocessing(cfg);