
%% SHPB Test Data Plots

clc
clear

time=zeros(125000,14);
v_incident=zeros(125000,14);
v_transmitted=zeros(125000,14);

for i=1:14
    filename=['Group4_Test',num2str(i),'.csv'];
    data=importdata(filename,',',16);
    time(:,i)=data.data(:,1);
    v_incident(:,i)=data.data(:,2);
    v_transmitted(:,i)=data.data(:,4);
end

for i=15:51
    filename=['data_',num2str(i),'.csv'];
    data=importdata(filename,',',16);
    time(:,i)=data.data(:,1);
    v_incident(:,i)=data.data(:,2);
    v_transmitted(:,i)=data.data(:,4);
end

color=[0 0 0; 1 0 0; 0 1 0; 0 0 1; .5 .5 .5];
    
% figure(1)
% hold on
%    
% for ii=1:5
%     plot(time(:,ii),v_incident(:,ii),'.-','Color',color(ii,:))
%     %plot(time(:,ii),v_transmitted(:,ii),'.','Color',color(ii,:))
% end
% legend({'Calibration','0.625 in Lead', '0.5 in Lead','0.375 in Lead','0.625 in Paper'})
% xlabel('Time (s)')
% ylabel('Voltage (v)')
% title('Calibration and Pulse Shaped Traces')
% set(gca,'FontSize',18)


%% Paper Plot

figure(1)
hold on
plot(time(:,1),v_incident(:,1),'-','Color',color(1,:))
plot(time(:,4),v_incident(:,4),'-','Color',color(4,:))

legend({'Calibration','Pulse Shaped'})
xlabel('Time [s]')
ylabel('Voltage [v]')
title('Calibration and Pulse Shaped Traces')
set(gca,'FontSize',14)
xlim([min(time(:,1)),max(time(:,1))])
%% Wave data analysis

speed = BarSpeed(5,2.4384);                                                %compute bar speed from recorded data
strength = zeros(1,46);                                                    %vector to store computed strength values
equilibrium = zeros(1,46);

for i=6:51                                                                 %each iteration of i will compute the strength value for a single dataset, set to just compute 1 to limit time
                                               
    dt = time(2,i)-time(1,i);                                              %time step in recorded data
    
%     figure(6)
%     plot(time(:,1),v_incident(:,i))
%     hold on
%     plot(time(:,i),v_transmitted(:,i))
    
    trigger = abs(mean(v_incident(1:100,i)));                              %average in signal noise of incident wave
    if trigger < 0.002
        trigger = 0.002;
    end
    
    count = 1;                                                             %used to index
    
    [a, b] = max(v_incident(:,i));                                         %a is max value of incident wave, b is index of this value
    [c, d] = min(v_incident(:,i));                                         %c is max value of reflected wave, d is index of this value

    index = floor(((1.2192*2)/speed)/dt);                                  %used bar speed to find beginning of reflected wave
        
    num = floor(((1.21920+.9652)/speed)/dt);                               %used bar speed to find index of transmitted wave rise
    
    for j=1:125000
        if v_incident(j,i) > 9.5*trigger                                   %average computed earlier used to trigger start of incident wave
            incident(count,1) = time(j,i);                                 %time values for incident wave
            incident(count,2) = v_incident(j,i);                           %voltage values for incident wave
            if count == 1
                mark = j;                                                  %index of trigger for incident wave
            end
            count = count+1;                                               %increased index
        end

    end
    
    n=1:1:length(incident);                                                %number of recorded data points for incident wave
    N = floor(length(incident)/2);                                         %Number of data points in important wave                                        
    T = incident(length(incident),1)-incident(1,1);                        %total elapsed time from start of incident wave to end
    rate = length(incident)/T;                                             %sampling rate of system
    k=1:1:N;
    w0 = (2*pi)/T;                                                         %natural frequency
    
    for m=1:length(n)
    
        reflect(m,1) = time(mark+m+index,i);                               %center index value found earlier used to compute start of reflected wave pulse, time values
        reflect(m,2) = v_incident(mark+m+index,i);                         %voltage values for reflected pulse
        transmitted(m,1) = time(mark+num+m,i);                             %center index value used to compute start of reflected pulse, time values
        transmitted(m,2) = v_transmitted(mark+num+m,i);                    %voltage values for reflected pulse
        
    end
    
    frequencies = zeros(1,N);                                              %all relevant frequency values
    Ck = zeros(1,N);                                   
    A=0.57106;
    B=0.43016;
    C=16.764;
    D=19.679;
    E=-5.8543;
    F=2.1532;
    radius = 0.01905/2;
    frequencies(1:end) = k(1:end).*(rate/N);                               %all relevant frequency values
    wavelength = speed./frequencies;                                       %wavelengths from frequencies
    
%     Ck(1:end) = speed.*(A+(B./((C*(radius./wavelength(1:end)).^4)+(D*(radius./wavelength(1:end)).^3)+(D*(radius./wavelength(1:end)).^2)+(F*(radius./wavelength(1:end)).^1.5)+1)));
    
    Ck = (k.*w0.*wavelength)/(2*pi);                                       %speed of each frequency
    
    incident_t = fft(incident(:,2))./(N);                                  %fft transform of incident wave
    reflect_t = fft(reflect(:,2))/(N/2);                                   %fft transform of reflected wave
    transmit_t = fft(transmitted(:,2))/(N/2);                              %fft transform of transmitted wave
    F_incident = zeros(1,length(n));
    F_reflect = zeros(1,length(n));
    F_transmit = zeros(1,length(n));
    for q=1:length(n)

        check = 0;
        for l=2:50
            phi(l) = k(l)*w0*(1.2192/Ck(l));
            check = check + real(incident_t(l))*cos((k(l-1)*w0*incident(q,1))-phi(l))+...  %forward/backward dispersion equation implementation for incident wave
                imag(incident_t(l))*sin((k(l-1)*w0*incident(q,1))-phi(l));
        end
        F_incident(q) = check+(real(incident_t(1))/2);
        
    end
    
    for x=1:length(n)
        summation = 0;
        for y=2:50
            phi(y) = k(y)*w0*(-0.9652/Ck(y));
            summation = summation + real(transmit_t(y))*cos((k(y-1)*w0*transmitted(x,1))-phi(y))+...  %forward/backward dispersion equation implementation of transmitted wave
                imag(transmit_t(y))*sin((k(y-1)*w0*transmitted(x,1))-phi(y));
        end
        F_transmit(x) = summation+(real(transmit_t(1))/2);
    end
    
    for x=1:length(n)
        summation = 0;
        for y=2:50
            phi(y) = k(y)*w0*(-1.2192/Ck(y));
            summation = summation + real(reflect_t(y))*cos((k(y-1)*w0*reflect(x,1))-phi(y))+...      %forward/backward dispersion equation implementation for reflected wave
                imag(reflect_t(y))*sin((k(y-1)*w0*reflect(x,1))-phi(y));
        end
        F_reflect(x) = summation+(real(reflect_t(1))/2);
    end
    
    t_incident = incident(:,1)+(1.2192/speed);                             %time index correction based on bar wave speed
    t_transmit = transmitted(:,1)-(0.9652/speed);
    t_reflect = reflect(:,1)-(1.2192/speed);
    
    [g,h]=max(F_transmit);
    e_transmit = abs(Volt2Strain(g));                                      %strains from max index of transmitted wave
    e_incident = (Volt2Strain(F_incident(h)));
    e_reflect = Volt2Strain(F_reflect(h));
    
    F1 = pi*((0.01905/2)^2)*70.3e9*e_incident+e_reflect;                   %Forces from strains
    F2 = pi*((0.01905/2)^2)*70.3e9*e_transmit;
    
    strength(i-5) = F2*(2/pi)*(1/(0.00635*0.01905));                       %strength values for each test
    equilibrium(i-5) = F2-F1;                                              %force equilibrium
    
end

%% Experiment 3 - Statistical Analysis

sigma = [strength(1:10) strength(12:30) strength(33:38)];
% excluded sample 11, 31, 32, and 39-46 due to sampling and incomplete wave
% forms
% Central Tendency:
sigma_bar = mean(sigma)  % mean
sigma_med = median(sigma)  % median

% Dispersion:
sigma_std = std(sigma)   % standard deviation
sigma_var = var(sigma)   % variance

% Weibull
sigma_dist = wblfit(sigma);
x_o = min(sigma)
b = sigma_dist(1)
m = sigma_dist(2)

figure(2)
wblplot(sigma)
xlabel('log(\sigma_t)')
title('Weibull Probability Dynamic Tensile Strength of Concrete')

%% Pressure vs Strength

P = [10.2, 10.2, 9, 9, 8, 8, 12.4, 12.4, 11.5, 11.9, 11.6, 12.4, 8.2, 8.1, 12, 12, 9.7, 8.4, 8.3, 11.1, 8.1, 8.2, 8.1, 8, 8.3, 8.3, 8.3, 8.3, 8.2, 8.4, 8.6, 8.6, 8.6, 8.6, 10.1];

p=polyfit(P,sigma/1E7,1);
Rs=corr(P,sigma,'type','Pearson');
fit=p(1)*P+p(2);

figure(3)
hold on
plot(P,sigma/1E7,'x')
plot(P,fit,'r','LineWidth',1.5)
xlabel('Gas Gun Pressure [psi]')
ylabel('\sigma_t [MPa]')
legend({'data',['Linear Fit']},'Location','NorthWest','FontSize',12)
title('Gas Gun Pressure vs Dyanmic Tensile Strength') 
set(gca,'FontSize',12)