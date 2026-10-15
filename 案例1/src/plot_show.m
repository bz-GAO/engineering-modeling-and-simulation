function plot_show(run_time, elite, final_best_cost, best_cost_history, avg_cost_history,volt, tempr, max_gen)
    % 结果统计与可视化模块
    fprintf('\n算法总耗时：%.2f 秒\n', run_time);
    fprintf('拟合的最低成本为: %.2f\n', final_best_cost);

    best_cost_index = find(elite == 1);
    best_cost_points = sum(elite);
    fprintf('所选取的拟合点有 %d 个，分别为 [%s]\n', best_cost_points, num2str(best_cost_index));

    %------- 绘图 ------
    close all;

    % 图1：标定成本图
    figure('Name', '标定成本图', 'Position', [100, 100, 500, 400]);
    plot(1:max_gen, best_cost_history, 'r-', 'LineWidth', 1.5);
    hold on;
    plot(1:max_gen, avg_cost_history, 'b--', 'LineWidth', 1.5);
    title('标定成本迭代曲线');
    xlabel('迭代代数');
    ylabel('最低标定成本');
    legend('每一代最小成本', '每一代平均成本');
    grid on;

    % 图2：标定函数对比图 (以样本1为例)
    figure('Name', '函数对比图', 'Position', [650, 100, 600, 400]);
    sample_idx = 1; % 提取第一个物理样本进行可视化展示
    true_V = volt(sample_idx, :);
    true_T = tempr(sample_idx, :);

    best_logic = logical(elite);
    pre_V = true_V(best_logic); % 提取测试点电压
    pre_T = true_T(best_logic); % 提取测试点温度

    % 获取插值生成的完整曲线
    fit_T = mycurvefitting(pre_V, pre_T, true_V);

    % 绘制真实的基准曲线与拟合曲线
    plot(true_V, true_T, 'k-', 'LineWidth', 1.5); hold on;
    plot(true_V, fit_T, 'r--', 'LineWidth', 1.5);
    % 散点标出到底选了哪几个测试点
    plot(pre_V, pre_T, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

    set(gca, 'XDir', 'reverse'); % 电压是随温度下降的，翻转X轴更符合直觉
    xlabel('输出电压 (mV)');
    ylabel('温度读数 (℃)');
    legend('真实特性曲线', '插值拟合曲线', '选定的标定测试点', 'Location', 'southeast');
    title(sprintf('样本 %d 上的插值效果对比 (选用 %d 个点)', sample_idx, best_cost_points));
    grid on;
end