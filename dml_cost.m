function J = dml_cost(theta_deg, Rhat, M, d, lambda)
% DML代价函数
% J(theta) = tr(P_A_perp * Rhat)
% P_A_perp = I - A*pinv(A)
%
% 输入：
%   theta_deg - 当前候选DOA集合（1xK 或 Kx1）
%   Rhat      - 样本协方差矩阵
%
% 输出：
%   J         - DML代价值，越小越好

    A = steering_ula(M, d, lambda, theta_deg);

    % 投影到A列空间上的投影矩阵
    P_A = A * pinv(A);

    % 正交补投影矩阵
    P_A_perp = eye(M) - P_A;

    % DML代价
    J = real(trace(P_A_perp * Rhat));
end