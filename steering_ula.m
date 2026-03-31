function A = steering_ula(M, d, lambda, theta_deg)
% 构造ULA导向矩阵
% 输入：
%   M         - 阵元数
%   d         - 阵元间距
%   lambda    - 波长
%   theta_deg - 角度（标量或向量，单位：度）
%
% 输出：
%   A         - M x K导向矩阵

    theta_deg = theta_deg(:).';          % 转成行向量
    m = (0:M-1).';                       % 阵元下标列向量

    % ULA导向矢量：a(theta) = exp(-j*2*pi*d/lambda * m * sin(theta))
    A = exp(-1i * 2*pi * d / lambda * m * sind(theta_deg));

    % 单位范数归一化：让每个导向矢量的2范数为1
    A = A / sqrt(M);
end