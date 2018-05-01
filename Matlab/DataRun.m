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
rate=[rate(1:3);zeros(length(StaticList),1)];
strength=[S_Ult(1:3,1);zeros(length(StaticList),1)];

for i=1:length(StaticList)
    [~,strength(i+length(1:3)),rate(i+length(1:3))]=QuasiStatic(StaticList{i});
end

K1=sqrt(pi*0.004)*strength/1E6;

figure
plot(log(rate),K1,'*')