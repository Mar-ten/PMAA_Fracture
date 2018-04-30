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

for i=1:length(List)
    [S_Ult(i,1),S_Ult(i,2)]=UltimateStrain(List{i});
end
