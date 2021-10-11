% FnAddMarkers
% read .Marker file and find the marker onsets
% Ehsan Eqlimi, WAVES group, Ghent University, November 2020
function [EEG,FinalEventName,FinalEventTimes]=FnAddMarkers (EEG,FilenameMarkers)
fid = fopen(FilenameMarkers,'rt') ;
Markers = textscan(fid,'%s','delimiter','\n') ;
%%%%% Find the indices of markers ( marker type and onset)
EventNameID=contains(Markers{1},'<Description>');
EventNames=Markers{1}(EventNameID);
EventTimes=Markers{1}(find(EventNameID)+1);
A = cellfun(@strfind,EventNames',repmat(cellstr('>'), 1,size(EventNames,1)),'UniformOutput',false);
B = cellfun(@strfind,EventNames',repmat(cellstr('<'), 1,size(EventNames,1)),'UniformOutput',false);
C = cellfun(@strfind,EventTimes',repmat(cellstr('>'), 1,size(EventTimes,1)),'UniformOutput',false);
D = cellfun(@strfind,EventTimes',repmat(cellstr('<'), 1,size(EventTimes,1)),'UniformOutput',false);
for y=1:length(EventNames)
    EEG.event(y).type=EventNames{y}(A{y}(1)+1:B{y}(2)-1);
    EEG.event(y).latency=str2num(EventTimes{y}(C{y}(1)+1:D{y}(2)-1));
    FinalEventName{y,1}=EventNames{y}(A{y}(1)+1:B{y}(2)-1);
    FinalEventTimes(y,1)=str2num(EventTimes{y}(C{y}(1)+1:D{y}(2)-1));
end
EEG.event=EEG.event(1:length(EventTimes));
EEG.urevent=EEG.event(1:length(EventTimes));