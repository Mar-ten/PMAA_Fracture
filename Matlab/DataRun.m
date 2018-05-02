List={'T0014All.csv';
    'T0015All.csv';
    'T0016All.csv';
    'T0017All.csv';
    'T0018All.csv';
    'T0019All.csv';
    'T0020All.csv';
    'T0021All.csv';
    'T0022All.csv'};

S_Ult=zeros(length(List),2);
rate=zeros(size(List));
for i=1:length(List)
    [S_Ult(i,1),S_Ult(i,2),rate(i)]=UltimateStrain(List{i});
end


StaticList={'Static0_5in-min.csv';'Static0_05in-min.csv';'Static0_005in-min.csv'};
rate=[rate;zeros(length(StaticList),1)];
strength=[S_Ult(:,1);zeros(length(StaticList),1)];

for i=1:length(StaticList)
    [~,strength(i+length(List)),rate(i+length(List))]=QuasiStatic(StaticList{i});
end

K1=sqrt(pi*0.004)*strength/1E6;
%% Strain Rate Plot


figure
h=plot(log(rate),K1,'d');
ax=gca;
grid on
xlabel('$log(\dot{\varepsilon}) \left[\frac{1}{s}\right]$','Interpreter','latex')
ylabel('$K_1 \left[MPa \sqrt{m}\right]$','Interpreter','latex')
h.MarkerSize=10;
ax.FontSize=16;
h.MarkerFaceColor=[1 0 1];
%ylim([0,8])
