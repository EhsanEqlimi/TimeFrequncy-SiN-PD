function OurLayout=FnEEGLayoutCreate(tmpeloc)
OurLayout=[];
for mm=1:length(tmpeloc)
OurLayout.pos(mm,1)=-tmpeloc(mm).Y;
OurLayout.pos(mm,2)=tmpeloc(mm).X;
OurLayout.width(mm,1)=0.2;
OurLayout.height(mm,1)=0.1;
OurLayout.label{mm}=tmpeloc(mm).labels;
end

OurLayout.pos(29,2)=-1.2;
OurLayout.pos(30,1)=-1.2;
OurLayout.pos(31,1)=+1.2;