%% demo_dml_doa_multi_source.m
% 多目标/多源 DML 测角完整示例
% 说明：
% 1) 阵列模型：远场、窄带、均匀线阵 ULA
% 2) 源信号模型：复高斯随机信号（也可以替换成确定性波形）
% 3) 噪声模型：独立同分布复高斯白噪声
% 4) DML准则：min tr(P_A_perp * Rhat)
% 5) 优化方式：Bartlett粗初始化 + DML交替投影/坐标下降细化

clear; clc; close all;

%% ========================= 1. 参数设置 =========================

% 阵列参数
M = 10;                          % 阵元数
c_light = 3*1e8;
fc = 77*1e9;
lambda = c_light/fc;             % 波长                    
d = lambda/2;                    % 阵元间距，常用lambda/2

% 信号参数
K = 3;                           % 源数（目标数）
theta_true = [-28, 27.7, 31.5];  % 真实DOA，单位：度
Nsnap = 500;                     % 快拍数
SNR_dB = 10;                     % 信噪比(dB)

% 搜索范围
search_range = [-60, 60];        % DML搜索角度范围（单位：度）

% 粗初始化参数（Bartlett）
init_grid_step = 0.5;            % 粗搜索步长
peak_min_sep = 3;                % 峰值最小间隔，避免多个初值挤在一起

% DML局部细化参数（多级分辨率）
refine_steps = [0.5, 0.1, 0.02]; % 每一级局部搜索的步长（单位：度）
refine_half_widths = [8, 2, 0.5];% 每一级搜索窗口半宽（单位：度）
max_iter_each_stage = 15;        % 每一级最多迭代次数

fprintf('================ DML多源测角示例 ================\n');
fprintf('阵元数 M = %d\n', M);
fprintf('源数   K = %d\n', K);
fprintf('真实DOA = %s (deg)\n', mat2str(theta_true));
fprintf('快拍数 N = %d\n', Nsnap);
fprintf('SNR = %.2f dB\n\n', SNR_dB);

%% ========================= 2. 生成多源信号 =========================
% 构造真实导向矩阵 A(theta_true)
A_true = steering_ula(M, d, lambda, theta_true);

% 生成多源复高斯随机信号（K x Nsnap）
% 每一行对应一个源，每一列对应一个快拍
S_true = (randn(K, Nsnap) + 1i*randn(K, Nsnap)) / sqrt(2);

% 这里将导向矢量做了单位范数归一化，所以每个源的阵列输出能量量级更稳定
% 噪声方差按SNR近似设置
sigma2_true = 10^(-SNR_dB/10);

% 生成复高斯白噪声（M x Nsnap）
Noise = sqrt(sigma2_true/2) * (randn(M, Nsnap) + 1i*randn(M, Nsnap));

% 生成阵列接收数据
% X = A*S + N
X = A_true * S_true + Noise;

%% ========================= 3. 样本协方差矩阵估计 =========================
Rhat = (X * X') / Nsnap;

% 如果要保证数值上严格Hermitian，可做一次对称化
Rhat = (Rhat + Rhat') / 2;

%% ========================= 4. 源数估计（可选） =========================
% 如果K未知，可以用MDL进行估计。可以默认使用已知K。
K_est = estimate_num_sources_mdl(Rhat, Nsnap, M);
fprintf('MDL估计得到的源数 K_est = %d\n', K_est);
K = K_est;

%% ========================= 5. Bartlett粗初始化 =========================
% 说明：
% DML本质是多维非线性优化，直接高维穷举代价很大。
% 所以先用Bartlett谱找K个粗峰值，作为DML初值，再做细化。

theta_grid = search_range(1) : init_grid_step : search_range(2);
P_bartlett = zeros(size(theta_grid));

for ii = 1:length(theta_grid)
    a = steering_ula(M, d, lambda, theta_grid(ii));
    % Bartlett空间谱（这里a已单位范数归一化）
    P_bartlett(ii) = real(a' * Rhat * a);
end

% 取K个峰值作为DML初始角度
theta_init = pick_topK_peaks(theta_grid, P_bartlett, K, peak_min_sep);
theta_init = sort(theta_init);

fprintf('Bartlett粗初始化 DOA = %s (deg)\n', mat2str(theta_init, 4));

%% ========================= 6. DML测角 =========================
% DML目标函数：
% J(theta) = tr(P_A_perp(theta) * Rhat)
% 其中 P_A_perp = I - A * pinv(A)
%
% 这里采用交替更新（坐标下降）：
% 固定其它K-1个角度，只更新其中1个角度；
% 轮流更新所有角度，多轮迭代直到收敛。

[theta_est, cost_history] = dml_ap_estimator( ...
    Rhat, K, M, d, lambda, theta_init, ...
    search_range, refine_steps, refine_half_widths, max_iter_each_stage);

theta_est = sort(theta_est);

%% ========================= 7. 回代估计信号与噪声方差 =========================
A_est = steering_ula(M, d, lambda, theta_est);

% 最小二乘估计源波形：S_hat = A^dagger * X
S_est = pinv(A_est) * X;

% 估计噪声方差：sigma2_hat = (1/M) * tr(P_A_perp * Rhat)
P_A_perp_est = eye(M) - A_est * pinv(A_est);
sigma2_est = real(trace(P_A_perp_est * Rhat)) / M;

%% ========================= 8. 输出结果 =========================
fprintf('\n================ 估计结果 ================\n');
fprintf('真实DOA     = %s (deg)\n', mat2str(sort(theta_true), 6));
fprintf('DML估计DOA  = %s (deg)\n', mat2str(theta_est, 6));
fprintf('噪声方差真值 = %.6f\n', sigma2_true);
fprintf('噪声方差估计 = %.6f\n', sigma2_est);

% 逐个输出误差（按排序后一一对应）
theta_true_sorted = sort(theta_true);
err = theta_est - theta_true_sorted;
fprintf('角度估计误差 = %s (deg)\n', mat2str(err, 6));

%% ========================= 9. 绘图 =========================
% 9.1 Bartlett谱（用于观察粗峰值）
P_bartlett_dB = 10*log10(abs(P_bartlett) / max(abs(P_bartlett)));

figure('Color', 'w');
plot(theta_grid, P_bartlett_dB, 'LineWidth', 1.5); grid on; hold on;
xlabel('Angle (deg)');
ylabel('Normalized Bartlett Spectrum (dB)');
title('Bartlett Spectrum (for initialization)');
ylim([-40, 0]);

for k = 1:K
    xline(theta_true(k), '--b', sprintf('True %.1f^\\circ', theta_true(k)), ...
        'LabelVerticalAlignment', 'middle');
end
for k = 1:K
    xline(theta_est(k), '--r', sprintf('DML %.2f^\\circ', theta_est(k)), ...
        'LabelVerticalAlignment', 'bottom');
end
legend('Bartlett Spectrum', 'Location', 'best');

% 9.2 DML代价函数下降过程
figure('Color', 'w');
plot(cost_history, '-o', 'LineWidth', 1.5, 'MarkerSize', 5); grid on;
xlabel('Iteration Index');
ylabel('DML Cost: tr(P_A^\perp Rhat)');
title('DML Cost Convergence');

% 9.3 给出最终阵列拟合残差能量
residual_energy = norm(X - A_est*S_est, 'fro')^2 / Nsnap;
fprintf('拟合残差平均能量 = %.6f\n', residual_energy);

%% ========================= 10. 可选：显示估计出的源波形能量 =========================
src_power_est = mean(abs(S_est).^2, 2);
fprintf('估计源信号平均功率 = %s\n', mat2str(src_power_est.', 6));









