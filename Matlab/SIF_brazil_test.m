filename='Static0_005in-min.csv';
data=importdata(filename,',',7);
Load=data.data(:,1);
Time=data.data(:,2);
Extension=data.data(:,3);
Stress = data.data(:,4);
Strain = data.data(:,5);





P = 100;
theta = 0;
[K_I,K_II,Mix_I,Mix_II] = SIF_Brazil(Load,theta);

plot(K_I,Load)