function [theta_est, cost_history] = dml_ap_estimator( ...
    Rhat, K, M, d, lambda, theta_init, ...
    search_range, refine_steps, refine_half_widths, max_iter_each_stage)
% 基于交替投影/坐标下降的DML求解器
%
% 输入：
%   Rhat                - 样本协方差矩阵
%   K                   - 源数
%   M, d, lambda        - 阵列参数
%   theta_init          - 初始DOA
%   search_range        - 总搜索范围 [min, max]
%   refine_steps        - 每一级局部搜索步长，例如 [0.5, 0.1, 0.02]
%   refine_half_widths  - 每一级局部搜索半宽，例如 [8, 2, 0.5]
%   max_iter_each_stage - 每一级最大迭代次数
%
% 输出：
%   theta_est           - 最终估计DOA
%   cost_history        - 代价收敛过程

    theta_est = sort(theta_init(:).');   % 当前DOA估计
    cost_history = [];

    if length(refine_steps) ~= length(refine_half_widths)
        error('refine_steps 和 refine_half_widths 的长度必须一致。');
    end

    for stage = 1:length(refine_steps)
        step = refine_steps(stage);
        half_width = refine_half_widths(stage);

        fprintf('\n---- DML 第 %d 级细化：步长 = %.4f deg, 窗口半宽 = %.4f deg ----\n', ...
            stage, step, half_width);

        for iter = 1:max_iter_each_stage
            theta_old = theta_est;

            % 轮流更新每一个角度
            for k = 1:K

                % 第1级的第1轮，可以直接用全局范围搜索，提高鲁棒性
                if stage == 1 && iter == 1
                    local_grid = search_range(1):step:search_range(2);
                else
                    left_bd  = max(search_range(1), theta_est(k) - half_width);
                    right_bd = min(search_range(2), theta_est(k) + half_width);
                    local_grid = left_bd:step:right_bd;
                end

                J_vals = zeros(size(local_grid));

                % 固定其它K-1个角度，仅扫描第k个角度
                for ii = 1:length(local_grid)
                    theta_try = theta_est;
                    theta_try(k) = local_grid(ii);

                    % 排序是为了保持角度顺序一致，方便后续输出
                    theta_try = sort(theta_try);

                    J_vals(ii) = dml_cost(theta_try, Rhat, M, d, lambda);
                end

                % 取使DML代价最小的角度
                [~, idx_best] = min(J_vals);
                theta_est(k) = local_grid(idx_best);

                % 每更新一个角度后重新排序
                theta_est = sort(theta_est);
            end

            % 记录当前代价
            curr_cost = dml_cost(theta_est, Rhat, M, d, lambda);
            cost_history(end+1) = curr_cost; 

            fprintf('Stage %d, Iter %2d: theta = %s, cost = %.8f\n', ...
                stage, iter, mat2str(theta_est, 6), curr_cost);

            % 收敛判据：本轮变化足够小
            if max(abs(theta_est - theta_old)) < step/5
                break;
            end
        end
    end
end