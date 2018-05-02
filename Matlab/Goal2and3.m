
angles=importdata('Angles2.xlsx');
beta=angles.data(:,2);
Theta=[angles.data(:,3),angles.data(:,4)];
SD=std(Theta,0,2);
Theta_Mean=angles.data(:,5);
K1 =zeros(size(beta)); K2=K1; N1=K1; N2=K1; 

a = 4/1000;
R=25/2/1000;
B=3/1000;

alpha=0:5:90;


for i= 1:length(beta)
    [K_I(i),K_II(i),N1(i),N2(i)]=SIF_Brazil(1,beta(i)); % Creates normalized values of K1 and K2 
end
loads =1.0e+03 *[2.3603;2.0797;1.7526;2.0330;2.0797;1.6124;1.7993;1.6124;1.7059];
% strength=strength(1:9);
% P=strength.*pi*B*R;
strength=loads/(pi*B*R);
for i= 1:length(beta)
    [K1(i),K2(i),N1(i),N2(i)]=SIF_Brazil(loads(i),beta(i)); % Creates normalized values of K1 and K2 
end
K_I=sqrt(pi*a).*cosd(beta).^2;
K_II=sqrt(pi*a).*cosd(beta).*sind(beta);

T_Max=-2*atand(K_I./(4.*K_II)-1/4*sqrt((K_I./K_II).^2+8));
T_Max(isnan(T_Max))=0;

KC=(K_I./4.*(3*cosd(T_Max./2)+cosd(3*T_Max./2)))-3/4*(sind(T_Max./2)+sind(3*T_Max./2)).*K_II;

% K1_Norm=K_I./KC;
% K2_Norm=K_II./KC;

% K1_Norm=K_I/K1_S;
% K2_Norm=K_II./K1_S;
% figure;
% plot(K2_Norm,K1_Norm)
% 
% K1_Norm=K_I./K1_S;
% K2_Norm=K_II./K2_S;


K1=strength.*sqrt(pi*a).*cosd(beta).^2;
K2=strength.*sqrt(pi*a).*cosd(beta).*sind(beta);


figure; hold on
plot(K_II/max(K_II),K_I/max(K_I),'-','LineWidth',2)
plot(K2/max(K2),K1/max(K1),'kd')


%plot(abs(K_II/max(K_II)),abs(K_I/max(K_I)),'-','LineWidth',2)
%plot(abs(K2/max(K2)),abs(K1/max(K1)),'kd')


xlabel('K_{II}/K_{II C}','FontSize',16)
ylabel('K_{I}/K_{I C}','FontSize',16)
legend({'Theoretical','Transmitted Bar'},'Location','best','FontSize',16)
set(gca,'FontSize',16)
grid on
%% Goal 3

figure; hold on
plot(beta,T_Max,'m','LineWidth',2)
h=errorbar(beta,Theta_Mean,SD,'kd');
ax=gca;
grid on
h.MarkerSize=10;
ax.FontSize=14;
legend({'Max Hoop Stress Criteria','Experimental Results'},'Location','best')
xlabel('Crack Orientation Angle (\beta)')
ylabel('Crack Kinking Angle (\theta")')