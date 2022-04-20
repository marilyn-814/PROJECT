% % 2维
% format long;
% Ndata = xlsread('./IGD结果/DTLZ6.xlsx','Sheet3','A3:A22'); % DTLZ3 NSGA2 3-obj
% Nave = xlsread('./IGD结果/DTLZ6.xlsx','Sheet3','B26'); % 对应的均值
% Nindex = find(Ndata>1);
% Ndata(Nindex,:) = Nave;
% 
% Mdata = xlsread('./IGD结果/DTLZ6.xlsx','Sheet3','C3:C22'); % DTLZ3 MOEA/D 3-obj
% Mave = xlsread('./IGD结果/DTLZ6.xlsx','Sheet3','B28'); % 对应的均值
% Mindex = find(Mdata>1);
% Mdata(Mindex,:) = Mave;
% 
% [h,p] = ttest(Ndata,Mdata);
% if p<0.05 % 显著性好
%     disp('YES');
% else
%     disp('NO'); % 表现相当
% end
% disp(p)


% 3维
format long;
Ndata = xlsread('./IGD结果/DTLZ7.xlsx','Sheet3','B3:B22'); % DTLZ3 NSGA2 3-obj
Nave = xlsread('./IGD结果/DTLZ7.xlsx','Sheet3','B27'); % 对应的均值
Nindex = find(Ndata>1);
Ndata(Nindex,:) = Nave;

Mdata = xlsread('./IGD结果/DTLZ7.xlsx','Sheet3','D3:D22'); % DTLZ3 MOEA/D 3-obj
Mave = xlsread('./IGD结果/DTLZ7.xlsx','Sheet3','B29'); % 对应的均值
Mindex = find(Mdata>1);
Mdata(Mindex,:) = Mave;

[h,p] = ttest(Ndata,Mdata);
if p<0.05 % 显著性好
    disp('YES');
else
    disp('NO'); % 表现相当
end
disp(p)
