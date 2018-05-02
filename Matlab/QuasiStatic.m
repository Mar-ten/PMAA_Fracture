function [ K1,S_Ult,rate ] = QuasiStatic( filename )

Data=importdata(filename,',',7); 
load=Data.data(:,1);
time=Data.data(:,2);
disp=Data.data(:,3);
stress=Data.data(:,4)*6894757.29; %convert Ksi to Pa
strain=Data.data(:,5);
a = 4/1000;

strainrate=diff(strain)./diff(time); 
rate=mode(strainrate); 

S_Ult=max(stress); 

K1=sqrt(pi*a)*S_Ult/1E6; 

end
