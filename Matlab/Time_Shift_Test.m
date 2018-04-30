% Time Shifting
close all
clear all

filename=['Group4_Test',num2str(8),'.csv'];
data=importdata(filename,',',16);
time=data.data(:,1);
v_incident=data.data(:,2);
v_transmitted=data.data(:,4);

% figure
% plot(time,v_incident,time,v_transmitted)

speed = 5.0732e+03; %Original Value
%speed = 5.05e+03;

L= 2.438/2;
L2=1.930/2;
dt=8.0000e-09;

time_shift=L/speed;
index_shift=round(time_shift/dt);

time_shift2=L2/speed;
index_shift2=round(time_shift2/dt);


%% Pad and Shift Data
pad=zeros(index_shift,1);
v_forward=[pad;v_incident;pad];  v_forward=TimeShift(v_forward,index_shift);   v_forward=downsample(v_forward,10);
v_back=[pad;v_incident;pad];     v_back=TimeShift(v_back,-index_shift);   v_back=downsample(v_back,10);
v_trans=[pad;v_transmitted;pad]; v_trans=TimeShift(v_trans,-index_shift2);   v_trans=downsample(v_trans,10);
v_diff=v_forward+v_back; 


% figure
% plot((1:length(v_forward))*dt,v_forward,(1:length(v_back))*dt,v_back,(1:length(v_back))*dt,v_diff,(1:length(v_back))*dt,v_trans)
% legend({'Incident','Reflected','Difference','Transmitted'},'FontSize',14);


%% Trim Data
Lcut=round(length(v_diff)/3);

v_forward(1:Lcut)=[]; v_forward(end-Lcut:end)=[];
v_back(1:Lcut)=[]; v_back(end-Lcut:end)=[];
v_diff(1:Lcut)=[]; v_diff(end-Lcut:end)=[];
v_trans(1:Lcut)=[]; v_trans(end-Lcut:end)=[];


%% Plots
figure;
plot((1:length(v_forward))*dt,v_forward,(1:length(v_back))*dt,v_back,(1:length(v_back))*dt,v_diff,(1:length(v_back))*dt,v_trans)
legend({'Incident','Reflected','Difference','Transmitted'},'FontSize',14);
xlim([0,max((1:length(v_forward))*dt)])
ylim([-0.3,0.3])
xlabel('Time [s]','FontSize',14)
ylabel('Voltage [v]','FontSize',14)
title('Aligned Waveforms')
set(gca,'FontSize',14)