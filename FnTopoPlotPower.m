function FnTopoPlotPower(MeanInducedPower,EEGChanLoc,MinVal,MaxVal,Type)
 figure,topoplot(MeanInducedPower, EEGChanLoc,'electrodes','labels','style','map','plotrad',.54)
 caxis([MinVal MaxVal]);
 h=colorbar;
 ylabel(h, Type)
 colormap('jet')

 
  