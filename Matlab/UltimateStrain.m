function [ F1max,F2max,rate ] = UltimateStrain( filename )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Time Shifting
% close all

data=importdata(filename,',',16);
time=data.data(:,1);
v_incident=data.data(:,2);
v_transmitted=data.data(:,4);
v_break=data.data(:,6);

speed = 5.0732e+03; %Original Value
%speed = 5.05e+03;
D=0.01905;          % Bar diameter = 3/4 inch
A=pi*((D/2)^2);     % cross sectional area of bar
E=70.3e9;           % elastic modulus of bar
t=0.00635;          %disk thickness
L1= 2.438/2;        % Length of incident bar
L2=1.930/2;         % Length of transmitted bar
dt=time(2)-time(1); % find the sample time

time_shift=L1/speed;
index_shift=round(time_shift/dt);

time_shift2=L2/speed;
index_shift2=round(time_shift2/dt);

%% Select Window
% 
% figure(1)
% plot(time,v_incident,time,v_transmitted,time,v_break);
% ax=gca;
% pts=getrect(ax);
% close(1)
% 
% T1=pts(1);
% DT=pts(3);
% i_1=round(T1/dt);
% i_2=round((T1+DT)/dt);
% 
% 
% v_incident(i_2:end)=[]; v_incident(1:i_1)=[];
% v_transmitted(i_2:end)=[]; v_transmitted(1:i_1)=[];
% v_break(i_2:end)=[]; v_break(1:i_1)=[];

%% Pad and Shift Data
pad=zeros(index_shift,1);
v_forward=[pad;v_incident;pad];  v_forward=TimeShift(v_forward,index_shift);   v_forward=downsample(v_forward,100);
v_back=[pad;v_incident;pad];     v_back=TimeShift(v_back,-index_shift);   v_back=downsample(v_back,100);
v_trans=[pad;v_transmitted;pad]; v_trans=TimeShift(v_trans,-index_shift2);   v_trans=downsample(v_trans,100);
v_break=[pad;v_break;pad]; v_break=downsample(v_break,100);
v_diff=v_forward+v_back;


% figure
% plot((1:length(v_forward))*dt,v_forward,(1:length(v_forward))*dt,v_back,(1:length(v_forward))*dt,v_diff,(1:length(v_forward))*dt,v_trans,(1:length(v_forward))*dt,v_break)
% legend({'Incident','Reflected','Difference','Transmitted','Conduction Gauge'},'FontSize',14);


%% Trim Data

Lcut=round(length(v_diff)/4);

% v_forward(1:Lcut)=[]; v_forward(end-Lcut:end)=[];
% v_back(1:Lcut)=[]; v_back(end-Lcut:end)=[];
% v_diff(1:Lcut)=[]; v_diff(end-Lcut:end)=[];
% v_trans(1:Lcut)=[]; v_trans(end-Lcut:end)=[];
% v_break(1:Lcut)=[]; v_break(end-Lcut:end)=[];


%% Plots
clipTime=(1:length(v_forward))*dt;

% figure;
% plot(clipTime,v_forward,clipTime,v_back,clipTime,v_diff,clipTime,v_trans,clipTime,v_break)
% legend({'Incident','Reflected','Difference','Transmitted','Conduction Gauge'},'FontSize',14);
% xlim([0,max((1:length(v_forward))*dt)])
% %ylim([-0.3,0.3])
% xlabel('Time [s]','FontSize',14)
% ylabel('Voltage [v]','FontSize',14)
% title('Aligned Waveforms')
% set(gca,'FontSize',14)

%% Strain

e_forward=Volt2Strain(v_forward);
e_back=Volt2Strain(v_back);
e_diff=Volt2Strain(v_diff);
e_trans=Volt2Strain(v_trans);

% figure;
% plot(clipTime,e_forward,clipTime,e_back,clipTime,e_diff,clipTime,e_trans)
% legend({'$\varepsilon_{forward}$','$\varepsilon_{back}$','$\varepsilon_{diff}$','$\varepsilon_{trans}$'},'Interpreter','latex','FontSize',14)
% xlabel('Time (s)')
% ylabel('Strain')
% title(filename)
%% Forces

F_1=A*E.*e_diff;
F_2=A*E.*e_trans;
F1max=max(F_1);
F2max=max(abs(F_2));
S_Ult1=2*max(F_1)/(pi*D*t);
S_Ult2=2*max(F_2)/(pi*D*t);

%% Get strain Rate

e_dot=abs(speed/0.0254*(-2*e_back));
rate=mean(e_dot);

% figure
% plot(clipTime,e_dot,clipTime,ones(size(clipTime))*rate)

end

