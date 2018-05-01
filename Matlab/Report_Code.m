


angles=importdata('Angles.xlsx');
beta=angles.data(:,2);
Theta=[angles.data(:,3),angles.data(:,4)];
SD=std(Theta,0,2);
Theta_Mean=angles.data(:,5);
K1=zeros(size(beta)); K2=K1; N1=K1; N2=K1; 

a = 4;
R=25/2;
B=3;

K_I=sqrt(pi*a).*cosd(beta).^2;
K_II=sqrt(pi*a).*cosd(beta).*sind(beta);

for i= 1:length(beta)
    [K1(i),K2(i),N1(i),N2(i)]=SIF_Brazil(1,beta(i)); % Creates normalized values of K1 and K2 
end


T_Max=2*atan2d(1,K_I./(4.*K_II)-1/4*sqrt((K_I./K_II).^2+8))-180;
T_Max(isnan(T_Max))=0;

K_I=sqrt(pi*a).*cosd(beta).^2;
K_II=sqrt(pi*a).*cosd(beta).*sind(beta);

KC=(1/4*(3*cosd(T_Max./2)+cosd(3*T_Max./2)).*K_I)-3/4*(sind(T_Max./2)+sind(3*T_Max./2)).*K_II;

K1_Norm=K_I./KC;
K2_Norm=K_II./KC;
figure; 
plot(K1_Norm,K2_Norm)
%% Goal 3


% K1_Norm=cosd(beta).^2./N1;
% K2_Norm=cosd(beta).*sind(beta)./N2;
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