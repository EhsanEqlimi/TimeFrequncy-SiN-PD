
function [EvokedPower,MeanEvokedPower]= FnEvokedPower(TFData,Method)
n=size(TFData,1);%trial number
% EvokedPower=[];
Xw=TFData;
TrialDim=1;
for Ch=1:size(Xw,2)
    ChanDat=squeeze(Xw(:,Ch,:,:));
    if strcmp(Method,'plv')
        EvokedPower(:,Ch,:)=abs(nanmean((ChanDat./abs(ChanDat)),TrialDim)).^2;
    elseif strcmp(Method,'itpc')
        EvokedPower(:,Ch,:)=abs(nanmean((ChanDat./abs(ChanDat)),TrialDim));
    elseif strcmp(Method,'itc')
        EvokedPower(:,Ch,:)=(abs(mean(ChanDat,TrialDim)).^2 )./(mean((abs(ChanDat).^2),TrialDim));
    elseif strcmp(Method,'itpcz')
        EvokedPower(:,Ch,:)= n*(abs(nanmean((ChanDat./abs(ChanDat)),TrialDim))).^2;
    end
end
MeanEvokedPower=squeeze(nanmean(EvokedPower,2));

