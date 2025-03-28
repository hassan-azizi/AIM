% marg_pos = readmatrix("CALF-20-CO2-case/CALF-20_margpos.csv");
% marg_pos_bin_edges = readmatrix("CALF-20-CO2-case/CALF-20_margpos_bin_edges.csv");
% bounds = readmatrix("CALF-20-CO2-case/CALF-20_params_bound.csv");
marg_pos = readmatrix("CALF-20_margpos.csv");
marg_pos_bin_edges = readmatrix("CALF-20_margpos_bin_edges.csv");
bounds = readmatrix("CALF-20_params_bound.csv");

lb = bounds(1, :);
ub = bounds(2, :);
np = size(marg_pos, 1);
num_a_params = 6;
param_names = ["a_{0}", "a_{1}", "a_{2}", "a_{3}", "a_{4}", "a_{5}",...
                "b_{0}", "b_{1}", "b_{2}", "b_{3}"];
% param_names = ["m_b", "b", "m_d", "d"];
param_units = ["K", "K(mol/kg)^{-1}", "K(mol/kg)^{-2}", "K(mol/kg)^{-3}", "K(mol/kg)^{-4}", "K(mol/kg)^{-5}",...
                "-", "K(mol/kg)^{-1}", "K(mol/kg)^{-2}", "K(mol/kg)^{-3}"];

leg_size = 12;
label_size = 18;
idx_scale = [];
scale_values = [10^4, 10^6];


theta_nominal =[-4.69E+03	2.63E+02	-8.55E+02	5.36E+02	-1.23E+02	9.441289...
                21.99285	0.688322	0.120564];

hist_values = marg_pos;
dist_params_norm = zeros(size(marg_pos, 1), 2);
dist_params_skew_norm = zeros(size(marg_pos, 1), 4);

for n=1:np-1
    dist_params_norm(n, :) = fit_norm(marg_pos_bin_edges(n, :), marg_pos(n, :), 2);
    dist_params_skew_norm(n, :) = fit_norm(marg_pos_bin_edges(n, :), marg_pos(n, :), 3);
end

for p=1:np-1
    fig1 = figure(p);
    fig1.Position = [100, 50, 500, 375];
    % subplot(3, 2, p);
    hold on
    % Plot marginal posterior histogram
    if ismember(p, idx_scale)
        idx = find(p==idx_scale);
        scale_value = scale_values(idx);
        h = histogram('BinEdges', marg_pos_bin_edges(p, :).*scale_value,...
              'BinCounts', marg_pos(p, :), 'Normalization', 'pdf');
        AIM_fit_line = xline(abs(theta_nominal(p)).*scale_value);
        xlabel(gca, sprintf("%s %c 10^{%d} (%s)", param_names(p), char(215),...
                int16(log10(scale_value)), param_units(p)), "FontSize",label_size);
        xlim(gca, [lb(p), ub(p)].*scale_value);
    else
        h = histogram('BinEdges', marg_pos_bin_edges(p, :),...
              'BinCounts', marg_pos(p, :), 'Normalization', 'pdf');
        x = get_centroid(marg_pos_bin_edges(p, :));
        % f = fun_norm(dist_params_norm(p, :), x);
        % f = skew_normal_pdf(dist_params_skew_norm(p, :), x);
        % plot(x, f, LineWidth=2.0);

        AIM_fit_line = xline(abs(theta_nominal(p)));
        xlabel(gca, sprintf("%s (%s)", param_names(p), param_units(p)),...
            "FontSize",label_size);
        xlim(gca, [lb(p), ub(p)]);
    end
    % h.FaceColor = "#0072BD";
    h.FaceColor = "#D95319";
    h.FaceAlpha = 0.7;
    AIM_fit_line.Layer = "top";
    % AIM_fit_line.Color = "#0072BD";
    AIM_fit_line.Color = "black";
    AIM_fit_line.LineStyle = "--";
    AIM_fit_line.LineWidth = 3.0;
    
    ax = gca;
    % xlabel(ax, sprintf("%s, (%s)", param_names(p), param_units(p)), "FontSize",label_size);
    ylabel(ax, sprintf("P(%s|D)", param_names(p)), "FontSize",label_size);
    set(ax, 'YTick', []);
    % xtickformat('%.1f')
    hold off
    % grid on
    ax.FontName = 'Deja Vu Sans';
    ax.FontSize = label_size;
    ax.GridAlpha = 0.3;
    box on;
    saveas(fig1, sprintf('CALF-20-CO2-case/%s.svg', param_names(p)))
end

% Q_ads UA
Q_ads_edges = marg_pos_bin_edges(1, :) .* 8.314e-03;
Q_ads_marg_pos = marg_pos(1, :) .* 8.314e-03;
Q_nominal = theta_nominal(1) .* 8.314e-03;
h = histogram('BinEdges', Q_ads_edges,...
              'BinCounts', Q_ads_marg_pos, 'Normalization', 'pdf');
AIM_fit_line = xline(abs(Q_nominal));
h.FaceColor = "#0072BD";
h.FaceAlpha = 0.7;
AIM_fit_line.Layer = "top";
% AIM_fit_line.Color = "#0072BD";
AIM_fit_line.Color = "black";
AIM_fit_line.LineStyle = "--";
AIM_fit_line.LineWidth = 3.0;
ax = gca;
low_lim = min([lb(1), ub(1)].*8.314e-03);
high_lim = max([lb(1), ub(1)].*8.314e-03);
xlim(ax, [low_lim, high_lim]);
xlabel(ax, sprintf("%s (%s)", "\Delta H_{ads}", "kJ/mol"));
ylabel(ax, sprintf("P(%s|D)", "\Delta H_{ads}"));
set(ax, 'YTick', []);
% set(ax, 'Fontsize', label_size, 'Linewidth', 1.2);
hold off
% grid on
ax.FontName = 'Deja Vu Sans';
ax.FontSize = label_size;
ax.GridAlpha = 0.3;
saveas(gcf, sprintf('CALF-20-CO2-case/Q_ads.svg'))
box on
close all;

%% Miscellenous Function
function fitted_params = fit_norm(bin_edges, bin_count, num_params)
    x = get_centroid(bin_edges);
    if num_params == 2
        lb = [min(x), 0]; % Lower bound (mu, sigma)
        ub = [max(x), 1]; % Upper bound (mu, sigma)
        fit_function = @(params) residual_norm(params, x, bin_count);
    else
        lb = [min(x), 0, -50]; % Lower bound (location, scale, shape)
        ub = [max(x), 1, -50]; % Upper bound (location, scale, shape)
        fit_function = @(params) residual_sk_norm(params, x, bin_count);
    end
    
    % Solver: lsqnonlin
    % problem = createOptimProblem('lsqnonlin', 'objective', @get_residual, 'xdata', x, 'ydata', y, 'lb', lb, 'ub', ub, 'x0', x0);
    % ms = MultiStart('Display', 'off', 'UseParallel', true);
    % rs = RandomStartPointSet('NumStartPoints', 2000);
    % sol = run(ms, problem, rs);
    
    options = optimoptions("ga", "Display", "off", "FunctionTolerance", 1e-14,...
                                "MaxTime", 5000, "MaxGenerations", 1000000);
    [fitted_params, residual, exitflag] = ga(fit_function, num_params, [], [], [], [], lb, ub,...
                 [], options);
    
    if num_params > 2
        delta = fitted_params(3)/sqrt(1+fitted_params(3)^2);
        mu_calc = fitted_params(1) + fitted_params(2)*delta*sqrt(2/pi());
        fitted_params = [fitted_params, mu_calc];
    end
    s = 0;
end

function x = get_centroid(bin_edges)
    % Calculate bin centroid locations
    x = zeros(1, length(bin_edges) - 1); % Bin centroid locations
    for i = 2:length(bin_edges)
	    x(i-1) = mean([bin_edges(i), bin_edges(i-1)]);
    end
end
function f = fun_norm(params, x)
    mu = params(1);
    sigma = params(2);
    SSE = (x-mu).^2;
    multiplier = 1/(sigma*sqrt(2*pi));
    exp_part = exp(-SSE/2/sigma^2);
    f = multiplier.*exp_part;
end

function [f] = skew_normal_pdf(params, x)
    % Unpack parameter values	
    e = params(1); % Location parameter
    w = params(2); % Scale parameter
    a = params(3); % Shape parameter
    
    % Caluclate vector of skew-normal PDF values (f)
    phi_pdf = (1/sqrt(2*pi))*exp(-((x-e).^2)/2/(w^2));
    phi_cdf = 0.5*(1 + erf(a*(x-e)/w/sqrt(2)));
    f = 2*phi_pdf.*phi_cdf/w;
end

function residual = residual_norm(params, x, y)
    f_pred = fun_norm(params, x);
    residual = sum(abs(f_pred - y).^2);
end

function residual = residual_sk_norm(params, x, y)
    f_pred = skew_normal_pdf(params, x);
    residual = sum(abs(f_pred - y).^2);
end

