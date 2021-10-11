
function MeanEvokedPower=FnFindEvokedPowerinSelTime(FreqRange,TimeRange,TFData,EvokedPower)

[~, Index1]=min(abs(repmat(TFData.freq,size(FreqRange,1))-FreqRange(:,1)),[],2);
[~, Index2]=min(abs(repmat(TFData.freq,size(FreqRange,1))-FreqRange(:,2)),[],2);
[~, Index3]=min(abs(repmat(TFData.time,size(TimeRange,1))-TimeRange(:,1)),[],2);
[~, Index4]=min(abs(repmat(TFData.time,size(TimeRange,1))-TimeRange(:,2)),[],2);
EvokedPower_Perm=permute(EvokedPower,[2 1 3]);
Power=EvokedPower_Perm(:,Index1:Index2,Index3:Index4 );
MeanEvokedPower=squeeze(nanmean(nanmean(Power,3),2));