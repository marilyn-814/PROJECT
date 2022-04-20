%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 筛选IGD最小时 对应的Pareto前沿，筛选合适的点 并画图
% MOEA/D 三维
% DTLZ1:在0.5内   DTLZ2-4：三个维度里都在1内   DTLZ7：x y 在1内；z 在6内
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj_fun = 'DTLZ7';
num_obj = 2;
% 理想Pareto前沿
pareto_front = Calculate_Pareto_Front(obj_fun, 10000, num_obj);
% 真实前沿
collected_result = MOEA_D      (obj_fun,num_obj);
allIGD = getIGD(collected_result);
minIGD = find(allIGD==min(allIGD)); %最小值对应的位置
minPF = collected_result(minIGD).non_dominated_front; %最小IGD对应的前沿
[x,y] = size(minPF);
% for i=1:x
%     for j=1:y
%         if y==3 && minPF(i,j)>1
%             minPF(i,j)=1;
%         elseif y==1 && minPF(i,j)>1
%             minPF(i,j)=1;
%         elseif y==2 && minPF(i,j)>1
%             minPF(i,j)=1;
%         end
%     end
% end
% for i=1:x
%     for j=1:y
%         if y==2 && minPF(i,j)>1
%             minPF(i,j)=1;
%         elseif y==1 && minPF(i,j)>1
%             minPF(i,j)=1;
%         end
%     end
% end

% scatter3(minPF(:,1),minPF(:,2),minPF(:,3),[],[0 0.4470 0.7410],'.');hold on;
% scatter3(pareto_front(:,1),pareto_front(:,2),pareto_front(:,3),[],[0.9290 0.6940 0.1250],'.');
% xlabel('f1');
% ylabel('f2');
% zlabel('f3');
scatter(minPF(:,1),minPF(:,2),[],[0 0.4470 0.7410],'o');hold on;
scatter(pareto_front(:,1),pareto_front(:,2),[],[0.9290 0.6940 0.1250],'.');

xlabel('f1');
ylabel('f2');
% legend('NSGA-II','PF');
% title(sprintf('MOEA/D on %d-objective %s \n',num_obj,obj_fun));
% title(sprintf('NSGA-II on %d-objective %s \n',num_obj,obj_fun));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 获取所有的IGD值
function allIGD = getIGD(collected_result)
    i = 1;
    allIGD = zeros(500,1);
    allIGD(1) = collected_result(1).IGD;
    while i < 500
        i = i+1;
        allIGD(i) = collected_result(i).IGD;
    end
end