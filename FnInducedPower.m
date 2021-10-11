function Power=FnInducedPower(TFData)
%TFData: Trial * Channel * freq * Time
Magnitude=abs(TFData);
Power=Magnitude.^2;
