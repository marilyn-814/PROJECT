function collected_result = NSGA2(obj_fun,num_obj)
    %%% collected_result 结构体 收集算法结果
    %%% collected_result(i).IGD -> IGD值 double
    %%% collected_result(i).non_dominated_front -> [] 对应的非支配前沿
    %%% obj_fun DLTZ1 2 3 4 5 6 7
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 非支配： 没有能够支配x的点，x称为非支配点  最大！
    % obj_fun = 'DTLZ2'; -> 测试函数
    % num_obj = 3;  -> 三维
    
    if num_obj == 2
        pop_size = 200;
    elseif num_obj == 3
        pop_size = 300;
    end
    
    num_vari = 10;
    evaluation = 0;
    max_gen = 500;
    generation = 1;
    lower_bound = zeros(1,num_vari);
    upper_bound = ones(1,num_vari);
    
    % t0 = cputime;

    % 获取理想的Pareto_front
    pareto_front = Calculate_Pareto_Front(obj_fun, 10000, num_obj);

    % 生成需要执行的200个坐标点  population
    pop_vari = repmat(lower_bound,pop_size,1) + rand(pop_size, num_vari).*repmat((upper_bound-lower_bound),pop_size,1); % 生成随机变量 200*10 矩阵
    pop_obj = feval(obj_fun,pop_vari,num_obj); % 生成对应的函数坐标 200*2

    % 获得 nondominated 的点
    % pop_size = 200
    % X 中的所有点和第一个点比较 小于第一个点标记为非支配
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m=size(pop_obj,1);
    Xmin=min(pop_obj);
    X1=pop_obj-Xmin(ones(m,1),:);     %make sure X1>=0;
    Xmean=mean(X1); 
    %sort X1 so that dominated points can be removed quickly
    [~,list]=sort(max(X1./(Xmean(ones(m,1),:)+max(Xmean)),[],2));
    Y=pop_obj(list,:);         
    membership = false(size(pop_obj,1),1);
    while numel(list) > 1
        k = list(1);
        X1 = Y(1,:); % pop_obj 的第一组坐标
        X_1 = repmat(X1,size(Y,1),1);
        X = Y - X_1; 
        nondominated = any(X<0, 2); % 非支配解在X中的序号
        membership(k) = all(any(X(nondominated,:)>0,2)); % 两个坐标都小于0
        Y = Y(nondominated,:);
        list = list(nondominated,:);
    end
    membership(list) = true;
    % 第一代的非支配前沿
    nondominatedf = pop_obj(membership,:);
    % 第一代的IGD
    IGD = mean(min(pdist2(pareto_front,nondominatedf),[],2));

    front_rank = Front_Rank(pop_obj,pop_size);
    crowd_distance = Crowding_Distance(pop_obj,front_rank);
    
    % 收集 collected_result 结果
    result.IGD = 0;
    result.non_dominated_front = [];
    collected_result = repmat(result,1,max_gen); % 初始化 结构体数组
    collected_result(generation).IGD = IGD;
    collected_result(generation).non_dominated_front = nondominatedf;
    evaluation = evaluation + pop_size;
    % 画图
    drawFig(num_obj,nondominatedf,pareto_front)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 遗传算法
    % parent selection using tournament selection
    while generation < max_gen
        k = 2; % 二元 tournament
        [~,rank] = sortrows([front_rank,-crowd_distance]);
        [~,rank] = sort(rank);
        parents  = randi(size(pop_vari,1),k,pop_size);
        [~,best] = min(rank(parents),[],1);
        index = parents(best+(0:pop_size-1)*k)';
        pop_parent = pop_vari(index,:);
        % crossover
        dis_c = 20;
        mu  = rand(pop_size/2,num_vari);
        parent1 = pop_parent(1:2:pop_size,:);
        parent2 = pop_parent(2:2:pop_size,:);
        beta = 1 + 2*min(min(parent1,parent2)-lower_bound,upper_bound-max(parent1,parent2))./max(abs(parent2-parent1),1E-6);
        alpha = 2 - beta.^(-dis_c-1);
        betaq = (alpha.*mu).^(1/(dis_c+1)).*(mu <= 1./alpha) + (1./(2-alpha.*mu)).^(1/(dis_c+1)).*(mu > 1./alpha);
        betaq = betaq.*(-1).^randi([0,1],pop_size/2,num_vari);
        offspring1 = 0.5*((1+betaq).*parent1 + (1-betaq).*parent2);
        offspring2 = 0.5*((1-betaq).*parent1 + (1+betaq).*parent2);
        pop_crossover = [offspring1;offspring2];
        % mutation
        dis_m = 20;
        pro_m = 1/num_vari;
        rand_var = rand(pop_size,num_vari);
        mu  = rand(pop_size,num_vari);
        deta = min(pop_crossover-lower_bound, upper_bound-pop_crossover)./(upper_bound-lower_bound);
        detaq = zeros(pop_size,num_vari);
        position1 = rand_var<=pro_m & mu<=0.5;
        position2 = rand_var<=pro_m & mu>0.5;
        detaq(position1) = ((2*mu(position1) + (1-2*mu(position1)).*(1-deta(position1)).^(dis_m+1)).^(1/(dis_m+1))-1);
        detaq(position2) = (1 - (2*(1-mu(position2))+2*(mu(position2)-0.5).*(1-deta(position2)).^(dis_m+1)).^(1/(dis_m+1)));
        pop_mutation = pop_crossover + detaq.*(upper_bound-lower_bound);
        pop_mutation  = max(min(pop_mutation,upper_bound),lower_bound);
        % calculate the objective values of offsprings
        % 更新后的种群目标
        pop_mutation_obj = feval(obj_fun,pop_mutation,num_obj);
        % environment selection
        pop_vari_inter = [pop_vari;pop_mutation];
        pop_obj_inter = [pop_obj;pop_mutation_obj];
        [front_rank,rank_num] = Front_Rank(pop_obj_inter,pop_size);
        crowd_distance = Crowding_Distance(pop_obj_inter,front_rank);
        next = front_rank < rank_num;
        last = find(front_rank==rank_num);
        [~,rank] = sort(crowd_distance(last),'descend');
        next(last(rank(1:pop_size-sum(next)))) = true;
        front_rank = front_rank(next);
        crowd_distance = crowd_distance(next);
        pop_vari = pop_vari_inter(next,:);
        pop_obj = pop_obj_inter(next,:);

        % 更新 evaluation number 和 generation number
        evaluation = evaluation + size(pop_mutation_obj,1);
        generation = generation + 1;

        % 新的非支配前沿
        non_dominated_front = Pareto_Set(pop_obj);
        % 获得 IGD
        IGD = mean(min(pdist2(pareto_front,non_dominated_front),[],2));
        
        collected_result(generation).IGD = IGD;
        collected_result(generation).non_dominated_front = non_dominated_front;
        
        drawFig(num_obj,non_dominated_front,pareto_front)

    end
    
%     CPU = cputime - t0;
%     disp(CPU);
    
end