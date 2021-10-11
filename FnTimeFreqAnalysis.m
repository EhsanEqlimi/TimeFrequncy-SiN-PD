function TF=FnTimeFreqAnalysis(Data,TFParam)
cfg=TFParam;
[TF]=ft_freqanalysis(cfg,Data);