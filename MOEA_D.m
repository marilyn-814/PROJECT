function collected_result = MOEA_D(obj_fun,num_obj)
    %%% collected_result 结构体 收集算法结果
    %%% collected_result(i).IGD -> IGD值 double
    %%% collected_result(i).non_dominated_front -> [] 对应的非支配前沿
    % obj_fun = 'DTLZ1';
    % num_obj = 3;
    
    max_gen = 500;
    num_vari = 10;
    evaluation = 0;
    generation = 1;
    num_neighbor = 20;
    lower_bound = 0*ones(1,num_vari);
    upper_bound = 1*ones(1,num_vari);
    
    t0 = cputime;
    
    % 生成均匀权重
    if num_obj == 2
        num_weight = 200;
        weight = UniformPoint(200,2); % m = 2 -> w1+w2 = 1
    elseif num_obj == 3
        % m = 3 -> w1+w2+w3 = 1
        num_weight = 300; % ->  200
        weight = UniformPoint(300,3);
    end   

    neighbor = pdist2(weight,weight);
     [~,neighbor] = sort(neighbor,2);
    neighbor = neighbor(:,1:num_neighbor);

    % the first population

    % pop_vari = lhsdesign(num_weight, num_vari,'criterion','maximin','iterations',1000).*(upper_bound - lower_bound) + lower_bound; % 100*10
    pop_vari = repmat(lower_bound,num_weight,1) + rand(num_weight, num_vari).*repmat((upper_bound-lower_bound),num_weight,1); % 生成随机变量 200*10 矩阵
    pop_obj = feval(obj_fun,pop_vari,num_obj); % 100*3
    ideal_point = min(pop_obj,[],1); % 1*3

    % 获取理想的Pareto_front
    pareto_front = Calculate_Pareto_Front(obj_fun, 10000, num_obj);
    % 第一代的IGD
    IGD = mean(min(pdist2(pareto_front,pop_obj),[],2));
    evaluation = evaluation + num_weight;
     % 收集 collected_result 结果
    result.IGD = 0;
    result.non_dominated_front = [];
    collected_result = repmat(result,1,max_gen); % 初始化 结构体数组
    collected_result(generation).IGD = IGD;
    collected_result(generation).non_dominated_front = pop_obj;
    
    % 画图
    % drawFig(num_obj,pop_obj,pareto_front,obj_fun,generation,evaluation,IGD)
    
    t0 = cputime;
    
    % 遗传算法
    while generation < max_gen
        for ii = 1:num_weight
            this_neighbor = neighbor(ii,:);
            parent = pop_vari(this_neighbor(randperm(num_neighbor,2)),:);
            % crossover
            dis_c = 20;
            mu  = rand(1,num_vari);
            parent1 = parent(1,:);
            parent2 = parent(2,:);
            beta = 1 + 2*min(min(parent1,parent2)-lower_bound,upper_bound-max(parent1,parent2))./max(abs(parent2-parent1),1E-6);
            alpha = 2 - beta.^(-dis_c-1);
            betaq = (alpha.*mu).^(1/(dis_c+1)).*(mu <= 1./alpha) + (1./(2-alpha.*mu)).^(1/(dis_c+1)).*(mu > 1./alpha);
            % the crossover is performed randomly on each variable
            betaq = betaq.*(-1).^randi([0,1],1,num_vari);
            offspring1 = 0.5*((1+betaq).*parent1 + (1-betaq).*parent2);
            offspring2 = 0.5*((1-betaq).*parent1 + (1+betaq).*parent2);
            crossover = [offspring1;offspring2];
            % mutation (ploynomial mutation)
            % dis_m is the distribution index of polynomial mutation
            dis_m = 20;
            pro_m = 1/num_vari;
            rand_var = rand(2,num_vari);
            mu  = rand(2,num_vari);
            deta = min(crossover-lower_bound, upper_bound-crossover)./(upper_bound-lower_bound);
            detaq = zeros(2,num_vari);
            position1 = rand_var<=pro_m & mu<=0.5;
            position2 = rand_var<=pro_m & mu>0.5;
            detaq(position1) = ((2*mu(position1) + (1-2*mu(position1)).*(1-deta(position1)).^(dis_m+1)).^(1/(dis_m+1))-1);
            detaq(position2) = (1 - (2*(1-mu(position2))+2*(mu(position2)-0.5).*(1-deta(position2)).^(dis_m+1)).^(1/(dis_m+1)));
            mutation = crossover + detaq.*(upper_bound-lower_bound); % repair 将所有的点控制在bound内
            mutation  = max(min(mutation,upper_bound),lower_bound);
            offspring = mutation(1,:);
            offspring_obj = feval(obj_fun,offspring,num_obj);
            % update the ideal point
            ideal_point =  min([pop_obj;offspring_obj],[],1);
            % update the populatuion
            g_old = max(abs(pop_obj(this_neighbor,:)-repmat(ideal_point,num_neighbor,1)).*weight(this_neighbor,:),[],2);
            g_new = max(repmat(abs(offspring_obj-ideal_point),num_neighbor,1).*weight(this_neighbor,:),[],2); % update 
            pop_vari(this_neighbor(g_new<=g_old),:) = repmat(offspring,sum(g_new<=g_old),1);
            pop_obj(this_neighbor(g_new<=g_old),:) = repmat(offspring_obj,sum(g_new<=g_old),1);
        end
        generation = generation+1;
        evaluation = evaluation+num_weight;
        
        IGD = mean(min(pdist2(pareto_front,pop_obj),[],2));
        collected_result(generation).IGD = IGD;
        collected_result(generation).non_dominated_front = pop_obj;
        % drawFig(num_obj,pop_obj,pareto_front,obj_fun,generation,evaluation,IGD);
        
    end
    
%     CPU = cputime - t0;
%     disp(CPU);
    
end