function K_est = estimate_num_sources_mdl(Rhat, Nsnap, M)
% MDL源数估计
% 适合在源数未知时使用
%
% 输入：
%   Rhat  - 样本协方差矩阵
%   Nsnap - 快拍数
%   M     - 阵元数
%
% 输出：
%   K_est - 估计的源数

    eigvals = sort(real(eig((Rhat + Rhat')/2)), 'descend');
    mdl = inf(1, M);

    for k = 0:(M-1)
        noise_eigs = eigvals(k+1:end);

        if isempty(noise_eigs)
            mdl(k+1) = inf;
            continue;
        end

        % 几何均值 / 算术均值
        g = exp(mean(log(noise_eigs + eps)));
        a = mean(noise_eigs);

        mdl(k+1) = -Nsnap * (M-k) * log(g/a) + 0.5 * k * (2*M-k) * log(Nsnap);
    end

    [~, idx] = min(mdl);
    K_est = idx - 1;
end