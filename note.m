%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 执行7个 DTLZ测试函数，收集不同的IDG结果 
% 图像
% 重复 20次 每次取最后一次得到的IGD 均值、方差
% 把每次运行的所有结果写入 DTLZn.xlsx 一个测试函数是一个文件
% sheet1 两个算法在二维和三维运行下的 每次运行的IGD值  generation NIGD_2 NIGD_3 MIGD_2 MIGD_3
% 1-500 genration 重复执行 20runs
% sheet2 记录最后以及最小的IGD对应的前沿坐标
% 2Nmin_x 2Nmin_y 3Nmin_x 3Nmin_y 3Nmin_z 2Mmin_x 2Mmin_y 3Mmin_x 3Mmin_y 3Mmin_z *500
% 2Nfin_x 2Nfin_y 3Nfin_x 3Nfin_y 3Nfin_z 2Mfin_x 2Mfin_y 3Mfin_x 3Mfin_y 3Mfin_z *500
% *20
% sheet3 通过最后迭代后最后一次的IGD 2N 3N 2M 3M *20 
%                         最小值
%                         均值
%                         方差
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj_fun = 'DTLZ5';
num_obj2 = 2; 
num_obj3 = 3;
max_iteration =  20;
iteration = 0;
generation = 500;

A_num = 500;
A = zeros(A_num,1);
for i=1:A_num
    A(i) = A(i) + i;
end

B_f = zeros(max_iteration,1);
C_f = zeros(max_iteration,1);
D_f = zeros(max_iteration,1);
E_f = zeros(max_iteration,1);

while iteration < max_iteration
    iteration = iteration + 1;
    B = NSGA2(obj_fun,num_obj2); % 2维NSGA2
    allB = getIGD(B);
    title1 = {'2Nfin_x','2Nfin_y'};
    title2 = {'2Nmin_x','2Nmin_y'};
    [finBT,minBT] = getTable2(B,generation,title1,title2);

    C = NSGA2(obj_fun,num_obj3); % 3维NSGA2
    allC = getIGD(C);
    title3 = {'3Nfin_x','3Nfin_y','3Nfin_z'};
    title4 = {'3Nmin_x','3Nmin_y','3Nmin_z'};
    [finCT,minCT] = getTable3(C,generation,title3,title4);

    D = MOEA_D(obj_fun,num_obj2); % 2维MOEA/D
    allD = getIGD(D);
    title5 = {'2Mfin_x','2Mfin_y'};
    title6 = {'2Mmin_x','2Mmin_y'};
    [finDT,minDT] = getTable2(D,generation,title5,title6);

    E = MOEA_D(obj_fun,num_obj3); % 3维MOEA/D
    allE = getIGD(E);
    title7 = {'3Mfin_x','3Mfin_y','3Mfin_z'};
    title8 = {'3Mmin_x','3Mmin_y','3Mmin_z'};
    [finET,minET] = getTable3(E,generation,title7,title8);

    % 写入sheet1 
    title = {'generation','NIGD_2','NIGD_3','MIGD_2', 'MIGD_3'}; % A B C D E
    data = [A allB allC allD allE]; 
    [m,p] = size(data);
    data_cell = mat2cell(data,ones(m,1),ones(p,1));%matrix转变成cell
    result = [title;data_cell];
    T = table(result);
    star_num = num2str(503*(iteration-1) + 1);
    end_num = num2str(503*(iteration-1) + 502);
    start = 'A';
    en = 'E';
    sym = ':';
    range = strcat(start,star_num,sym,en,end_num);
    writetable(T,'DTLZ1.xlsx','Sheet',1,'Range',range);%保存文件
    
    % 写入sheet2
    % AB  CDE  FG  HIJ  KL  MNO  PQ  RST
    % finBT finCT finDT finET || minBT minCT minDT minET
    start1 = 'A';
    en1 = 'B';
    start2 = 'C';
    en2 = 'E';
    start3 = 'F';
    en3 = 'G';
    start4 = 'H';
    en4 = 'J';
    start5 = 'K';
    en5 = 'L';
    start6 = 'M';
    en6 = 'O';
    start7 = 'P';
    en7 = 'Q';
    start8 = 'R';
    en8 = 'T';
    sym = ':';
    % star_num = num2str(503*(iteration-1) + 1);
    % end_num = num2str(503*(iteration-1) + 502);
    range1 = strcat(start1,star_num,sym,en1,end_num);
    range2 = strcat(start2,star_num,sym,en2,end_num);
    range3 = strcat(start3,star_num,sym,en3,end_num);
    range4 = strcat(start4,star_num,sym,en4,end_num);
    range5 = strcat(start5,star_num,sym,en5,end_num);
    range6 = strcat(start6,star_num,sym,en6,end_num);
    range7 = strcat(start7,star_num,sym,en7,end_num);
    range8 = strcat(start8,star_num,sym,en8,end_num);
    % finBT finCT finDT finET || minBT minCT minDT minET
    writetable(finBT,'DTLZ1.xlsx','Sheet',2,'Range',range1);%保存文件
    writetable(finCT,'DTLZ1.xlsx','Sheet',2,'Range',range2);%保存文件
    writetable(finDT,'DTLZ1.xlsx','Sheet',2,'Range',range3);%保存文件
    writetable(finET,'DTLZ1.xlsx','Sheet',2,'Range',range4);%保存文件
    writetable(minBT,'DTLZ1.xlsx','Sheet',2,'Range',range5);%保存文件
    writetable(minCT,'DTLZ1.xlsx','Sheet',2,'Range',range6);%保存文件
    writetable(minDT,'DTLZ1.xlsx','Sheet',2,'Range',range7);%保存文件
    writetable(minET,'DTLZ1.xlsx','Sheet',2,'Range',range8);%保存文件
    
    % 写入sheet3
    finalB = B(generation).IGD;
    B_f(iteration) = finalB;
    finalC = C(generation).IGD;
    C_f(iteration) = finalC;
    finalD = D(generation).IGD;
    D_f(iteration) = finalD;
    finalE = E(generation).IGD;
    E_f(iteration) = finalE;
end

% 写入sheet3
% A B C D 1-22
% 2N 3N 2M 3M  
title_all = {'NSGA22','NSGA23','MOEAD2','MOEAD3'};
data_all = [B_f C_f D_f E_f];
[m_all,p_all] = size(data_all);
cell_all = mat2cell(data_all,ones(m_all,1),ones(p_all,1));%matrix转变成cell
resultall = [title_all;cell_all];
T_all = table(resultall); 

% 最小值
N2_min = min(B_f);
N3_min = min(C_f);
M2_min = min(D_f);
M3_min = min(E_f);
% 筛选IGD值 去掉所有大于0.1的值
indB = find(B_f<0.1);
B_f = B_f(indB);
numB = size(B_f,1);

indC = find(C_f<0.1);
C_f = C_f(indC);
numC = size(C_f,1);

indD = find(D_f<0.1);
D_f = D_f(indD);
numD = size(D_f,1);

indE = find(E_f<0.1);
E_f = E_f(indE);
numE = size(E_f,1);

% 均值
N2_ave = sum(B_f)/numB;
N3_ave = sum(C_f)/numC;
M2_ave = sum(D_f)/numD;
M3_ave = sum(E_f)/numE;
% 方差
N2_dev = 0;N3_dev = 0;M2_dev = 0;M3_dev = 0;
for i=1:numB
    N2_dev = (N2_dev*(i-1)+(B_f(i)-N2_ave)^2)/i;
end
for i=1:numC
    N3_dev = (N3_dev*(i-1)+(C_f(i)-N3_ave)^2)/i;
end
for i=1:numD
    M2_dev = (M2_dev*(i-1)+(D_f(i)-M2_ave)^2)/i;
end
for i=1:numE
    M3_dev = (M3_dev*(i-1)+(E_f(i)-M3_ave)^2)/i;
end

title_ana = {'min','Average','Variance'};
data_ana = [N2_min N2_ave N2_dev;N3_min N3_ave N3_dev;M2_min M2_ave M2_dev;M3_min M3_ave M3_dev];
[m_ana,p_ana] = size(data_ana);
cell_ana = mat2cell(data_ana,ones(m_ana,1),ones(p_ana,1));
result_ana = [title_ana;cell_ana];
T_ana = table(result_ana); 
writetable(T_all,'DTLZ1.xlsx','Sheet',3);
writetable(T_ana,'DTLZ1.xlsx','Sheet',3,'Range','A24:C29');



function allIGD = getIGD(collected_result)
    i = 1;
    allIGD = zeros(500,1);
    allIGD(1) = collected_result(1).IGD;
    while i < 500
        i = i+1;
        allIGD(i) = collected_result(i).IGD;
    end
end


function [Tfin,Tmin] = getTable2(B,generation,title1,title2)
    allB = getIGD(B);
    [minB,indexB] = min(allB); 
    % 2Nfin_x 2Nfin_y 最后的
    Nfin2 = B(generation).non_dominated_front;
    Nfinx_2 = Nfin2(:,1);
    Nfiny_2 = Nfin2(:,2);
    Nfin2d = [Nfinx_2 Nfiny_2];
    % title1 = {'2Nfin_x','2Nfin_y'};
    [m1,p1] = size(Nfin2d);
    data_cell1 = mat2cell(Nfin2d,ones(m1,1),ones(p1,1));%matrix转变成cell
    result1 = [title1;data_cell1];
    Tfin = table(result1); 
    % 2Nmin_x 2Nmin_y 最小的
    Nmin2 = B(indexB).non_dominated_front;
    Nminx_2 = Nmin2(:,1);
    Nminy_2 = Nmin2(:,2);
    Nmin2d = [Nminx_2 Nminy_2];
    % title2 = {'2Nmin_x','2Nmin_y'};
    [m1,p1] = size(Nmin2d);
    data_cell2 = mat2cell(Nmin2d,ones(m1,1),ones(p1,1));%matrix转变成cell
    result2 = [title2;data_cell2];
    Tmin = table(result2); 
end

function [Tfin,Tmin] = getTable3(B,generation,title1,title2)
    allB = getIGD(B);
    [minB,indexB] = min(allB); 
    % 2Nfin_x 2Nfin_y 最后的
    Nfin3 = B(generation).non_dominated_front;
    Nfinx_3 = Nfin3(:,1);
    Nfiny_3 = Nfin3(:,2);
    Nfinz_3 = Nfin3(:,3);
    Nfin3d = [Nfinx_3 Nfiny_3 Nfinz_3];
    [m1,p1] = size(Nfin3d);
    data_cell1 = mat2cell(Nfin3d,ones(m1,1),ones(p1,1));%matrix转变成cell
    result1 = [title1;data_cell1];
    Tfin = table(result1); 
    % 2Nmin_x 2Nmin_y 最小的
    Nmin3 = B(indexB).non_dominated_front;
    Nminx_3 = Nmin3(:,1);
    Nminy_3 = Nmin3(:,2);
    Nminz_3 = Nmin3(:,3);
    Nmin3d = [Nminx_3 Nminy_3 Nminz_3];
    [m1,p1] = size(Nmin3d);
    data_cell2 = mat2cell(Nmin3d,ones(m1,1),ones(p1,1));%matrix转变成cell
    result2 = [title2;data_cell2];
    Tmin = table(result2); 
end
