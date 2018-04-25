%% SHPB Test Data Plots

clc
clear
close all

time=zeros(125000,14);
v_incident=zeros(125000,14);
v_transmitted=zeros(125000,14);

for i=1:5
    filename=['Group4_Test',num2str(i),'.csv'];
    data=importdata(filename,',',16);
    time(:,i)=data.data(:,1);
    v_incident(:,i)=data.data(:,2);
    v_transmitted(:,i)=data.data(:,4);
end

% for i=15:51
%     filename=['data_',num2str(i),'.csv'];
%     data=importdata(filename,',',16);
%     time(:,i)=data.data(:,1);
%     v_incident(:,i)=data.data(:,2);
%     v_transmitted(:,i)=data.data(:,4);
% end

color=[0 0 0; 1 0 0; 0 1 0; 0 0 1; .5 .5 .5];
    
figure(1)
hold on
   
for ii=1:5
    A = downsample(time(:,ii),100);
    B=downsample(v_incident(:,ii),100);
    %x=smooth(A);
    %y = smooth(B);
    plot(A,B,'-','Color',color(ii,:))
    %plot(time(:,ii),v_transmitted(:,ii),'.','Color',color(ii,:))
end
legend({'Calibration','0.625 in Lead', '0.5 in Lead','0.375 in Lead','0.625 in Paper'})
title('Pulse Shaper and Calbration wave forms for SHPB')
xlabel('Time [s]');
ylabel('Volts [mV]')
    