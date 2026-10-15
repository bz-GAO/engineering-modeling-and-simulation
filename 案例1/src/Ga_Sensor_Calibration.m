%---主程序---

clear; clc; close all; tic;
rng(2026); % 固定随机种子，便于复现实验

script_dir = fileparts(mfilename('fullpath'));
data_file = fullfile(script_dir, '..', 'data', 'dataform_train2026.csv');
if ~isfile(data_file)
    error('缺少训练数据：%s。请参照 data/README.md 放置课程数据。', data_file);
end

%---常量定义与变量的初始化---
pop_size = 20;     % 已改参数
Q = 70;            
max_gen = 10;     
pc = 0.9;           
pm_start = 0.1;    
pm_end = 0.02;      

fprintf('\n开始读入数据...\n');
minput = readmatrix(data_file);
[M, gene_len] = size(minput);
volt = minput(2:2:end, :);
tempr = minput(1:2:end, :); % 奇数位为温度，偶数位为电压

best_cost_history = zeros(max_gen, 1);
avg_cost_history = zeros(max_gen, 1);
%---种群初始化---
%一定程度上加快运行速度，在早期运行得到一定低成本后，仅在4-12个点的范围内考虑
pop = zeros(pop_size, gene_len);
for i = 1:pop_size
    num_ones = randi([4, 12]);
    idx = randperm(gene_len, num_ones);
    pop(i, idx) = 1;
end

fprintf('开始遗传算法迭代...\n');

%---主循环---
for gen = 1 : max_gen
    costs = zeros(pop_size, 1);
    
    %--- 成本评估 ---
    for i = 1:pop_size
        costs(i) = eval_cost(pop(i, :), volt, tempr, Q);
    end
    
    %--- 适应度计算 ---
    fitness = zeros(pop_size, 1);
    valid_mask = costs < 10000; 
    
    if any(valid_mask)
        max_valid = max(costs(valid_mask));
        % 仅对合法个体进行评估，并使用平方级放大，拉开优劣差距
        % 采用最大值减去成本，使得小成本的具有大适应度，同时拉开相对差距
        fitness(valid_mask) = (max_valid - costs(valid_mask) + 10).^2; 
        fitness(~valid_mask) = 1e-6; % 非法个体给予极小概率
    else
        fitness(:) = 1; % 兜底保护
    end

    %---精英保留之取出---
    [min_cost, best_index] = min(costs);
    best_cost_history(gen) = min_cost;
    elite = pop(best_index, :);
    
    %--- 计算并记录平均成本 ---
    if any(valid_mask)
        avg_cost_history(gen) = mean(costs(valid_mask)); % 核心：只对"合法"个体求平均
    else
        avg_cost_history(gen) = mean(costs); % 极端情况兜底
    end
   
    
    elite = pop(best_index, :);

    %--- 轮盘赌选择 ---
    prob = fitness / sum(fitness);
    cumsum_prob = cumsum(prob);
    new_pop = zeros(pop_size, gene_len);
    for i = 1:pop_size
        r = rand();
        idx = find(cumsum_prob >= r, 1, 'first');
        if isempty(idx), idx = pop_size; end
        new_pop(i, :) = pop(idx, :);
    end
    pop = new_pop;


    %--- 交叉 --- 
    for i = 1:2:pop_size-1
        if rand() < pc
            cross_point = randi([1, gene_len-1]);
            temp = pop(i+1, cross_point+1:end);
            pop(i+1, cross_point+1:end) = pop(i, cross_point+1:end);
            pop(i, cross_point+1:end) = temp;
        end
    end

    %--- 变异 (混合策略：前期降温突变，后期升温滑移) ---
    phase_split = max_gen * 0.5; % 设定前50%为大范围搜索期
    
    if gen < phase_split
        %单点降温突变
        pm_early_start = 0.1;
        pm_early_end = 0.02;
        current_pm = pm_start - (gen / phase_split) * (pm_start - pm_end);
        
        for i = 1:pop_size
            for j = 1:gene_len
                if rand() < current_pm
                    pop(i, j) = 1 - pop(i, j);
                end
            end
        end
    else
        % 基因升温滑移
        pm_slide_start = 0.15;
        pm_slide_end = 0.3;
        progress = (gen - phase_split) / (max_gen - phase_split);
        current_slide_pm = pm_slide_start + progress * (pm_slide_end - pm_slide_start);
        
        for i = 1:pop_size
            climb_points = find(pop(i, :) == 1);
            num_points = length(climb_points);
            
            % 跳板机制：仅允许接近最优解的个体（6或7个点）进行滑移试探
            if (num_points == 6 || num_points == 7) && rand() < current_slide_pm
                pt_idx = randi([1, num_points]);
                slide_idx = climb_points(pt_idx);
                
                if rand() < 0.5
                    target_idx = slide_idx - 1;
                else
                    target_idx = slide_idx + 1;
                end
                
                if target_idx >= 1 && target_idx <= gene_len && pop(i, target_idx) == 0
                    pop(i, slide_idx) = 0;
                    pop(i, target_idx) = 1;
                end
            else
                % 兜底保护：非6或7点的个体保留 0.01 的常规突变率，避免过早陷入死水
                for j = 1:gene_len
                    if rand() < 0.01
                        pop(i, j) = 1 - pop(i, j);
                    end
                end
            end
        end
    end
    
    %--- 柔性边界保护与精英存入 ---
    % 在早期强制要求左右两端(各5个点)必须至少有1个点，保证测定有效,在后期则取消，防止卡死
   if gen < max_gen * 0.8
        for i = 1:pop_size
            if sum(pop(i, 1:5)) == 0, pop(i, randi([1, 5])) = 1; end
            if sum(pop(i, 86:90)) == 0, pop(i, randi([86, 90])) = 1; end
        end
   end
    
    pop(1, :) = elite;
end

%--- 爬山算法的最后优化  ---
[elite, final_best_cost] = hill_climbing_search(elite, best_cost_history(max_gen), volt, tempr, Q, gene_len);

%--- 结果统计与绘图  ---
run_time = toc;
plot_show(run_time, elite, final_best_cost, best_cost_history, avg_cost_history,volt, tempr, max_gen);
