function tmpeloc=FnEEGChanLocCreate (Elec)

% % Elec(20,:)=[];
% % Elec(32,:)=[];
% % Elec.Var1{32}='32';
for chan=1:size(Elec,1)
EEGChanLoc(chan).type='EEG';
EEGChanLoc(chan).labels=Elec.Var2{chan};
EEGChanLoc(chan).label=Elec.Var2{chan};
EEGChanLoc(chan).sph_phi_besa=(Elec.Var4(chan));
EEGChanLoc(chan).sph_theta_besa=(Elec.Var3(chan));
end
writelocs( EEGChanLoc, 'EEGChanLoc.elp','filetype','besa' );
[tmpeloc labels Th Rd indices] = readlocs( 'EEGChanLoc.elp','filetype' ,'besa');
