function MeanInducedPower=FnFindInducedPowerinSelTime(FreqRange,TimeRange,TFData,MeanPowerTrials)

[~, Index1]=min(abs(repmat(TFData.freq,size(FreqRange,1))-FreqRange(:,1)),[],2);
[~, Index2]=min(abs(repmat(TFData.freq,size(FreqRange,1))-FreqRange(:,2)),[],2);
[~, Index3]=min(abs(repmat(TFData.time,size(TimeRange,1))-TimeRange(:,1)),[],2);
[~, Index4]=min(abs(repmat(TFData.time,size(TimeRange,1))-TimeRange(:,2)),[],2);
PowerPreSent=MeanPowerTrials(:,Index1:Index2,Index3:Index4 );
MeanInducedPower=squeeze(nanmean(nanmean(PowerPreSent,3),2));