function theta_init = pick_topK_peaks(theta_grid, spectrum, K, min_sep_deg)
% 从Bartlett谱中选择K个峰值作为初始化
%
% 输入：
%   theta_grid  - 扫描角度网格
%   spectrum    - 对应谱值
%   K           - 需要选取的峰数
%   min_sep_deg - 峰值最小间隔
%
% 输出：
%   theta_init  - 初始角度估计

    spectrum = spectrum(:).';
    theta_grid = theta_grid(:).';

    % 找局部峰值
    is_peak = false(size(spectrum));
    for i = 2:length(spectrum)-1
        if spectrum(i) >= spectrum(i-1) && spectrum(i) >= spectrum(i+1)
            is_peak(i) = true;
        end
    end

    peak_idx = find(is_peak);

    % 如果局部峰值太少，则退化为全局排序
    if isempty(peak_idx)
        [~, sort_idx] = sort(spectrum, 'descend');
        peak_idx = sort_idx;
    else
        [~, order] = sort(spectrum(peak_idx), 'descend');
        peak_idx = peak_idx(order);
    end

    theta_init = [];

    % 先按最小间隔挑选峰值
    for idx = peak_idx
        cand = theta_grid(idx);
        if isempty(theta_init) || all(abs(cand - theta_init) >= min_sep_deg)
            theta_init(end+1) = cand; %#ok<AGROW>
        end
        if length(theta_init) >= K
            break;
        end
    end

    % 如果还不够K个，用全局最大值补齐
    if length(theta_init) < K
        [~, sort_idx] = sort(spectrum, 'descend');
        for idx = sort_idx
            cand = theta_grid(idx);
            if isempty(theta_init) || all(abs(cand - theta_init) >= min_sep_deg/2)
                theta_init(end+1) = cand; %#ok<AGROW>
            end
            if length(theta_init) >= K
                break;
            end
        end
    end

    theta_init = sort(theta_init(1:K));
end