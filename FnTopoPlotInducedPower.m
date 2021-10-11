function FnTopoPlotPower (MeanInducedPower,EEGChanLoc,MinVal,MaxVal)
 figure,topoplot(MeanInducedPower, EEGChanLoc,'electrodes','labels','style','map','plotrad',.54)
 caxis([MinVal MaxVal]);
 h=colorbar;
 ylabel(h, 'Induced pre-sentence alpha power')
 colormap('jet')

 
  