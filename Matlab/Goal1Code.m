close all

List={'T0014All.csv';
    'T0015All.csv';
    'T0016All.csv';
    'T0017All.csv';
    'T0018All.csv';
    'T0019All.csv';
    'T0020All.csv';
    'T0021All.csv';
    'T0022All.csv'};


a = 4/1000;
R=25/2/1000;
B=3/1000;

angles=importdata('Angles2.xlsx');
beta=angles.data(:,2);
beta=[beta;0;0;0];
F_Ult=zeros(length(List),2);
rate=zeros(size(List));
K1=zeros(length(List),1);
for i=1:length(List)
    [F_Ult(i,1),F_Ult(i,2),rate(i)]=UltimateStrain(List{i});
end


StaticList={'Static0_5in-min.csv';'Static0_05in-min.csv';'Static0_005in-min.csv'};
rate=[rate;zeros(length(StaticList),1)];

%load=zeros(size(rate));


load=[F_Ult(:,2);zeros(length(StaticList),1)]; % choose 1 for differences signal, 2 for transmitted signal

%load(1:9)=strength(1:9)*pi*R*B;


for i=1:length(StaticList)
    [~,strength(i+length(List)),rate(i+length(List)),load(i+length(List))]=QuasiStatic(StaticList{i});
end


% P=strength*pi*R*B;

for i= 1:length(beta)
    [K1(i),K2(i),N1(i),N2(i)]=SIF_Brazil(load(i),beta(i)); % Creates normalized values of K1 and K2 
end

%K1=sqrt(pi*0.004)*strength/1E6;
%% Strain Rate Plot 1
figure
rate2=[rate(1:3);rate(end-2:end)];  % Takes values from tests with angle of 0
K12=[K1(1:3);K1(end-2:end)]./1E6;    % Takes values from tests with angle of 0
h=plot(log10(rate2),K12,'d');
ax=gca;
grid on
xlabel('$log(\dot{\varepsilon}) \left[\frac{1}{s}\right]$','Interpreter','latex')
ylabel('$K_1 \left[MPa \sqrt{m}\right]$','Interpreter','latex')
h.MarkerSize=10;
ax.FontSize=16;
h.MarkerFaceColor=[1 0 1];
%ylim([0,8])
