function drawFig(num_obj,pop_obj,pareto_front)
    %%% num_obj ; 算法获取的前沿 ; 理想前沿 ; obj_fun ; generation ; evaluation ; IGD
    %%% 画图 理想的Pareto_front 和 算法计算得到的 Pareto_front
    if num_obj == 2
        scatter(pop_obj(:,1),pop_obj(:,2),[],[0 0.4470 0.7410],'o');hold on;
        scatter(pareto_front(:,1),pareto_front(:,2),[],[0.9290 0.6940 0.1250],'.');
        xlabel('f1');
        ylabel('f2');
%         scatter(pop_obj(:,1),pop_obj(:,2),'ro'); hold on;
%         scatter(pareto_front(:,1),pareto_front(:,2),'b.');
        % title(sprintf('NSGA-II on %d-objective %s \n generation: %d, evaluations: %d, IGD: %0.4g',num_obj,obj_fun,generation,evaluation,IGD));
        % title(sprintf('MOEA/D on %d-objective %s \n generation: %d, evaluations: %d, IGD: %0.4g',num_obj,obj_fun,generation,evaluation,IGD));
        drawnow;hold off;
    elseif num_obj == 3
        scatter3(pop_obj(:,1),pop_obj(:,2),pop_obj(:,3),[],[0 0.4470 0.7410],'o');hold on;
        scatter3(pareto_front(:,1),pareto_front(:,2),pareto_front(:,3),[],[0.9290 0.6940 0.1250],'.');
        xlabel('f1');
        ylabel('f2');
        zlabel('f3');
%         scatter3(pop_obj(:,1),pop_obj(:,2),pop_obj(:,3),'ro'); hold on;
%         scatter3(pareto_front(:,1),pareto_front(:,2),pareto_front(:,3),'b.');
        % title(sprintf('NSGA-II on %d-objective %s \n generation: %d, evaluations: %d, IGD: %0.4g',num_obj,obj_fun,generation,evaluation,IGD));
        % title(sprintf('MOEA/D on %d-objective %s \n generation: %d, evaluations: %d, IGD: %0.4g',num_obj,obj_fun,generation,evaluation,IGD));
        view(135,30);drawnow;hold off;
    end
end