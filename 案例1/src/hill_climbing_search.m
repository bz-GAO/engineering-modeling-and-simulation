function [best_gene, best_cost] = hill_climbing_search(elite, current_best_cost, volt, tempr, Q, gene_len)
    % 爬山算法优化模块
    fprintf('调用爬山算法优化中...');
    best_gene = elite;
    best_cost = current_best_cost;
    need_improve = true;
    climb_iter = 0;

    while need_improve
        need_improve = false;
        climb_iter = climb_iter + 1;

        climb_points = find(best_gene == 1);

        for k = 1 : length(climb_points)
            pt = climb_points(k); % 获取第k个选定点的索引
            
            % 左查找
            if pt > 1 && best_gene(pt - 1) == 0 % 左有空位
                temp_gene = best_gene;
                temp_gene(pt) = 0;
                temp_gene(pt - 1) = 1;
                temp_cost = eval_cost(temp_gene, volt, tempr, Q);

                if temp_cost < best_cost
                    best_cost = temp_cost;
                    best_gene = temp_gene;
                    need_improve = true;
                    break; % 有更好解就跳出当下for，重新筛选
                end         
            end

            % 右查找
            if pt < gene_len && best_gene(pt + 1) == 0 % 右有空位
                temp_gene = best_gene;
                temp_gene(pt) = 0;
                temp_gene(pt + 1) = 1;
                temp_cost = eval_cost(temp_gene, volt, tempr, Q);

                if temp_cost < best_cost
                    best_cost = temp_cost;
                    best_gene = temp_gene;
                    need_improve = true;
                    break;
                end
            end
        end
    end
end