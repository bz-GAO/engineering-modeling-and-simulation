% =========================================================================
% 多节点声纳系统中同步时钟机制的可靠性评估和系统优化
% 基于变步长事件驱动与能力分桶算法的蒙特卡罗仿真
% =========================================================================

clear; clc; close all;
rng(2026); % 固定随机种子，便于比较不同节点规模

%% 1. 全局参数设置
k = 6;                  % 系统正常工作所需的最少有效节点数
S = 500000;             % 蒙特卡罗样本数 (推荐 10^5 级别以保证平滑度)
w = 25000;              % 目标可靠性考核时间 (hours)
N_range = 6:20;         % 冗余节点规模的遍历范围

mean_Tf_results = zeros(size(N_range));%储存每一次的结果，最终取最优
reliability_results = zeros(size(N_range));

%% 2. 核心映射矩阵 (O(1) 差分更新)
% 定义节点状态枚举值: 1=PF, 2=SO, 3=DM, 4=MO, 5=DN, 6=FB
%PF:完全正常，SO:只能从节点，DM:失效或者主节点，MO:仅主节点，DN:失效，FB:阻塞总线
% 映射矩阵 STATE_MAP: 行对应切换器 A 的状态 (1=A0, 2=A1, 3=A2, 4=A3)
%                     列对应切换器 B 的状态 (1=B0, 2=B1, 3=B2)
%具体做法是按照A、B的状态从STATE_MAP里查找，时间复杂度很低，做变步长变化很方便
STATE_MAP = [
    1, 4, 2;   % A0 + [B0, B1, B2] -> [PF, MO, SO]
    2, 6, 2;   % A1 + [B0, B1, B2] -> [SO, FB, SO]
    3, 4, 5;   % A2 + [B0, B1, B2] -> [DM, MO, DN]
    5, 5, 5    % A3 + [B0, B1, B2] -> [DN, DN, DN]
];

%% 3. 外层循环：遍历节点规模 n，做不同规模下的比较
tic;
for n_idx = 1:length(N_range) %对于可能的节点数都进行试验
    n = N_range(n_idx);
    T_f_array = zeros(1, S); % 记录每个样本的首次失效时间
    
    % 内层
    for s = 1:S
        T_f = 250000; % 把最大寿命作为默认值，防永生比较方便
        
        % -----------------------------------------------------------
        % 步骤 A: 预先生成所有元器件的寿命和故障模式
        % -----------------------------------------------------------
        % 生成发生时间 (指数分布均值 1/lambda)
        % 逆变换抽样，避免依赖 Statistics and Machine Learning Toolbox
        time_A = -55000 * log(rand(1, n));
        time_B = -220000 * log(rand(1, n));
        life_comp = [time_A, time_B]; %存放故障时间，轮取最小值后更新的inf可考虑用250000代替
        
        % 预先判定发生故障时的突变模式 (A:1,2,3  B:1,2)
        r_A = rand(1, n);
        mode_A = ones(1, n);                  % A1: 0.6，此处先默认生成1，后面再按概率覆盖2和3
        mode_A(r_A > 0.6 & r_A <= 0.75) = 2;  % A2: 0.15
        mode_A(r_A > 0.75) = 3;               % A3: 0.25
        
        r_B = rand(1, n);
        mode_B = ones(1, n);                  % B1: 0.65
        mode_B(r_B > 0.65) = 2;               % B2: 0.35
        
        mode_all = [mode_A, mode_B];         
        %存放故障类型，根据life_comp取出的最小时间和索引范围决定最新的变化
        
        % -----------------------------------------------------------
        % 步骤 B: 初始化状态追踪器与分桶
        % -----------------------------------------------------------
        % 0表示健康。后续加上偏移量1用于匹配 STATE_MAP 的行列索引(State_map中1对应正常索引)
        state_A = zeros(1, n); 
        state_B = zeros(1, n);
        
        % Buckets 索引: 1=PF, 2=SO, 3=DM, 4=MO, 5=DN, 6=FB
        buckets = zeros(1, 6);
        buckets(1) = n; % 初始时，所有 n 个正常节点都在 PF 桶中
        
        % -----------------------------------------------------------
        % 步骤 C: 事件驱动循环
        % -----------------------------------------------------------
        while true
            % 变步长找下一个故障
            [current_time, ev_idx] = min(life_comp);
                        
            if current_time >= 250000 % 显然当处理到我们二次更新的250000时也会退出
                break; %懒得单独写while变量了，直接做个break判断
            end
            
            % 标记该事件已被处理
            life_comp(ev_idx) = 250000; 
            
            % 判断是哪个节点的哪个切换器坏了，更新
            if ev_idx <= n
                node_idx = ev_idx;
                old_state_val = STATE_MAP(state_A(node_idx)+1, state_B(node_idx)+1);%后面要旧组里面删去，所以这个旧值也得取出
                state_A(node_idx) = mode_all(ev_idx); % 更新A
                new_state_val = STATE_MAP(state_A(node_idx)+1, state_B(node_idx)+1);
            else
                node_idx = ev_idx - n;
                old_state_val = STATE_MAP(state_A(node_idx)+1, state_B(node_idx)+1);
                state_B(node_idx) = mode_all(ev_idx); % 更新B
                new_state_val = STATE_MAP(state_A(node_idx)+1, state_B(node_idx)+1);
            end
            
            % 更新能力分桶
            buckets(old_state_val) = buckets(old_state_val) - 1;
            buckets(new_state_val) = buckets(new_state_val) + 1;
            
            % 提取各种点数量 (此处不关系具体序号考虑数量，5_DN直接忽略)
            Q_PF = buckets(1); Q_SO = buckets(2); Q_DM = buckets(3); 
            Q_MO = buckets(4); Q_FB = buckets(6);
            
            % -------------------------------------------------------
            % 步骤 D: 判断
            % -------------------------------------------------------
            system_failed = false;
            
            % 绝对数量不足
            % double((Q_MO + Q_DM) > 0)
            % 代表只要有主节点备选，就至多只可能有一个节点贡献，这里考虑极端的点数情况，直接考虑这个会被用
            % 分别对应总线阻塞、主节点冲突(MO不能静默)、没有主节点备选、有效节点数不足k
            if Q_FB >= 1 || Q_MO >= 2 || (Q_PF + Q_MO + Q_DM) == 0 || (Q_PF + Q_SO + double((Q_MO + Q_DM) > 0) < k)
                system_failed = true;
                

            % k-1从节点和一些主节点，但如果其中有PF被选成主，显然会爆
            elseif Q_MO == 0 && Q_PF >= 1 && Q_DM >= 1 && (Q_PF + Q_SO == k - 1)
                P_pick_DM = Q_DM / (Q_DM + Q_PF);% 按比例抽取
                if rand() >= P_pick_DM
                    % 运气不佳，系统挑中了 PF 当主节点，导致从节点名额 -1，跌破 k 阈值
                    system_failed = true;
                end
            end
            
            % 结算寿命
            if system_failed
                T_f = current_time;
                break;
            end
            
        end % 结束单个样本的 while 事件流
        
        T_f_array(s) = T_f;
    end % 结束 S 次大循环
    
    % 记录当前节点规模 n 的统计指标
    mean_Tf_results(n_idx) = mean(T_f_array);
    reliability_results(n_idx) = sum(T_f_array >= w) / S;
    
    fprintf('规模 n=%2d 仿真完成 | 平均寿命: %8.2f h | R(25000): %.4f\n', ...
            n, mean_Tf_results(n_idx), reliability_results(n_idx));
end
toc;

% --- 寻找最优的节点数 n ---
[max_MTTF, idx_MTTF] = max(mean_Tf_results);
optimal_n_MTTF = N_range(idx_MTTF);

[max_R, idx_R] = max(reliability_results);
optimal_n_R = N_range(idx_R);

fprintf('\n=== 最优结果分析 ===\n');
fprintf('使得平均寿命最大的节点数: n = %d, E[T_f] = %.2f\n', optimal_n_MTTF, max_MTTF);
fprintf('使得系统可靠性最大的节点数: n = %d, R(25000) = %.4f\n', optimal_n_R, max_R);

%% 4. 可视化图表生成
figure('Position', [100, 100, 900, 400]);

% 绘制平均工作寿命曲线
subplot(1, 2, 1);
plot(N_range, mean_Tf_results, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
grid on;
title('节点数 n 与系统平均工作寿命 E[T_f] 的关系');
xlabel('节点数 n');
ylabel('平均工作寿命 E[T_f] (hours)');

% 绘制可靠性曲线
subplot(1, 2, 2);
plot(N_range, reliability_results, '-s', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'Color', 'r');
grid on;
title('节点数 n 与系统可靠性 R(25000) 的关系');
xlabel('节点数 n');
ylabel('系统可靠性 R(25000)');
