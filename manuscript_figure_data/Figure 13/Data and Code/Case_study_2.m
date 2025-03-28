clc;
clear;
close all;
%% Loading Data and Consistency
data_1 = readmatrix("./Data/Isotherm_273.csv");
data_2 = readmatrix("./Data/Isotherm_298.csv");
data_3 = readmatrix("./Data/Isotherm_323.csv");

temperature_data = [273, 298, 323];

new_data = data_consistency([], data_1(:, 1:2));
new_data = data_consistency(new_data, data_2(:, 1:2));
new_data = data_consistency(new_data, data_3(:, 1:2));

pressure_data = new_data(:, 1:2:end);
loading_data = new_data(:, 2:2:end);

p_lin = rmmissing(reshape(pressure_data, [], 1));
q_lin = rmmissing(reshape(loading_data, [], 1));
N = size(p_lin, 1);
P_mean = mean(log(p_lin));
% P_mean = mean(log(p_lin));
%
%% Isotherm Function
iso_model = 'Virial';
iso_struc = isotherm_fit_opt(iso_model, loading_data, pressure_data);
iso_fun = iso_struc.fun;
%
%% Fitted Parameters
% fit_params = [-4.845673e+03, 9.654406e02, -1.223853e03, 5.838586e02, -1.221194e2, 9.281350,...
%                  22.501361, -1.547888, 1.277603, -0.149109];
fit_params = [-4.686E+03	2.63E+02	-8.55E+02	5.36E+02	-1.23E+02	9.441289...
                21.99285	0.688322	0.120564];
M = size(fit_params, 2);    % Number of params in DS-Langmuir Function
num_a = 6;
%
%% RMSE Calculation
% RMSE_values = zeros(length(method_array), 1);
% r2_values = zeros(length(method_array), 1);
% q_pred = zeros(length(method_array), N);
% 
% for i=1:length(method_array)
%     y = iso_fun(fit_params(i, :), P);
%     SSE = sum((q-y).^2);
% 
%     RMSE_values(i) = sqrt(SSE/(N-M));
%     r2_values(i) = 1 - SSE/(sum((q-q_mean).^2));
%     q_pred(i, :) = y;
% end

y = iso_fun(fit_params, loading_data, temperature_data, num_a);
% y = exp(y);
y_lin = rmmissing(reshape(y, [], 1));

SSE = sum((log(p_lin)-y_lin).^2);
% SSE = sum((pressure_data-y).^2);

RMSE_values = sqrt(SSE/(N-M));
r2_values = 1 - SSE/(sum((log(p_lin)-P_mean).^2));
%
%% Plots
fig1 = figure(1);
fig1.Position = [100, 50, 850, 750];
legend_str = {};
hold on
for i=1:length(temperature_data)
    s_p = scatter(rmmissing(pressure_data(:, i)), rmmissing(loading_data(:, i)), 'filled');
    s_p.SizeData = 100;
    s_p.MarkerEdgeColor = 'k';
    s_p.LineWidth = 0.5;
    s_p.MarkerFaceAlpha = 0.60; 
    % s_p.MarkerFaceColor = 'k'; 
    s_p.LineWidth = 1.0;
    legend_str = [legend_str, {sprintf('GCMC (%.0fK)', temperature_data(i))}];
    % legend_str = [legend_str, {sprintf('%.0fK', temperature_data(i))}];
end

for i=1:length(temperature_data)
    l_p = plot(rmmissing(exp(y(:, i))), rmmissing(loading_data(:, i)));
    l_p.LineWidth = 2.5;
    % p(i).Color = color_array(i);
    legend_str = [legend_str, {sprintf('HeatFit (%.0fK)', temperature_data(i))}];
end
hold off

ax1 = gca();
xlabel(ax1, "Pressure (Pa)", FontSize=26);

ylabel(ax1, "CO_2 Uptake (mol/kg)", FontSize=26)
legend(ax1, legend_str, FontSize=20, Location="best", NumColumns=1);

grid on
box on
consistent_color(ax1);
ax1.FontName = 'Deja Vu Sans';
ax1.GridAlpha = 0.3;
ax1.FontSize = 22;
ax1.XTick = [10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6];
% ax1.XLim = [0, 600];
% ax1.YLim = [0, 2.0];
xscale(ax1, "log");
% yscale(ax1, "log");
ytickformat('%.1f')
xtickangle(ax1, 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig2 = figure(2);
fig2.Position = [100, 50, 850, 750];
legend_str = {};
hold on
for i=1:length(temperature_data)
    l_p = plot(rmmissing(loading_data(:, i)), rmmissing(y(:, i)));
    l_p.LineWidth = 2.5;
    % p(i).Color = color_array(i);

    s_p = scatter(rmmissing(loading_data(:, i)), rmmissing(log(pressure_data(:, i))),  'filled');
    s_p.SizeData = 100;
    % s_p.MarkerEdgeColor = 'k';
    s_p.LineWidth = 0.6;
    s_p.MarkerFaceAlpha = 0.60; 
    % s_p.MarkerFaceColor = 'k'; 
    s_p.LineWidth = 1.0;
    legend_str = [legend_str, {sprintf('%.0fK', temperature_data(i)), ''}];
end
hold off

ax1 = gca();
xlabel(ax1, "Gas Uptake (mol/kg)", FontSize=26);

ylabel(ax1, "ln P (Pa)", FontSize=26)
legend(ax1, legend_str, FontSize=20, Location="best", NumColumns=1);

grid on
box on
consistent_color(ax1);
ax1.FontName = 'Deja Vu Sans';
ax1.GridAlpha = 0.3;
ax1.FontSize = 22;
% ax1.XLim = [0, 600];
ax1.YLim = [0, 14];
% xscale(ax1, "log");
ytickformat('%.1f')
xtickangle(ax1, 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GCMC heat
data_1 = readmatrix("./Data/heat_273.csv");
data_2 = readmatrix("./Data/heat_298.csv");
data_3 = readmatrix("./Data/heat_323.csv");

Q_ads = 8.314 .* polyval(flip(fit_params(1:5)), q_lin)./1e03;  % kJ/mol
fig3 = figure(3);
fig3.Position = [100, 50, 850, 750];
sz = 20;
hold on
% e1 = errorbar(data_1(:, 1), -abs(data_1(:, 2)), data_1(:, 3));
% e2 = errorbar(data_2(:, 1), -abs(data_2(:, 2)), data_2(:, 3));
% e3 = errorbar(data_3(:, 1), -abs(data_3(:, 2)), data_3(:, 3));

e1 = errorbar(data_1(:, 1), abs(data_1(:, 2)), data_1(:, 3));
e2 = errorbar(data_2(:, 1), abs(data_2(:, 2)), data_2(:, 3));
e3 = errorbar(data_3(:, 1), abs(data_3(:, 2)), data_3(:, 3));

e1.CapSize = sz;
e1.MarkerFaceColor = "#EDB120";
e1.Color = "#EDB120";
e1.LineWidth = 2.0;
e1.DisplayName = "273K";

e2.CapSize = sz;
e2.MarkerFaceColor = "#D95319";
e2.Color = "#D95319";
e2.LineWidth = 2.0;
e2.DisplayName = "298K";

e3.CapSize = sz;
e3.MarkerFaceColor = "#0072BD";
e3.Color = "#0072BD";
e3.LineWidth = 2.0;
e3.DisplayName = "323K";

% s_p = scatter(data_1(:, 1), -abs(data_1(:, 2)), 'filled');
% s_p.SizeData = 110;
% s_p.MarkerEdgeColor = 'k';
% s_p.LineWidth = 0.4;
% s_p.MarkerFaceAlpha = 0.80; 
% s_p.MarkerFaceColor = '#EDB120'; 
% s_p.LineWidth = 1.0;  
% 
% s_p = scatter(data_2(:, 1), -abs(data_2(:, 2)), 'filled');
% s_p.SizeData = 110;
% s_p.MarkerEdgeColor = 'k';
% s_p.LineWidth = 0.4;
% s_p.MarkerFaceAlpha = 0.80; 
% s_p.MarkerFaceColor = '#D95319'; 
% s_p.LineWidth = 1.0;  
% 
% s_p = scatter(data_3(:, 1), -abs(data_3(:, 2)), 'filled');
% s_p.SizeData = 110;
% s_p.MarkerEdgeColor = 'k';
% s_p.LineWidth = 0.4;
% s_p.MarkerFaceAlpha = 0.80; 
% s_p.MarkerFaceColor = '#0072BD'; 
% s_p.LineWidth = 1.0;  


s_p = scatter(q_lin, abs(Q_ads),  'filled');
s_p.SizeData = 120;
s_p.MarkerEdgeColor = 'k';
s_p.LineWidth = 0.5;
s_p.MarkerFaceAlpha = 0.90; 
s_p.MarkerFaceColor = 'r'; 
s_p.LineWidth = 1.0;   
% plot(q_lin, Q_ads,"LineWidth", 2.5, "Color", "r");
hold off

ax1 = gca();
xlabel(ax1, "CO_2 Uptake (mol/kg)", FontSize=26);

% ylabel(ax1, "Heat of Adsorption (kJ/mol)", FontSize=26)
ylabel(ax1, strcat(char(8722), '\DeltaH_{ads} (kJ/mol)'), 'Interpreter', 'tex', FontSize=26)
% legend(ax1, {'GCMC-273K', 'GCMC-298K','GCMC-323K', 'HeatFit'}, FontSize=20, Location="best", NumColumns=1);
legend(ax1, {'GCMC (273K)', 'GCMC (298K)','GCMC (323K)', 'HeatFit'}, FontSize=20, Location="best", NumColumns=1);

grid on
box on
ax1.FontName = 'Deja Vu Sans';
ax1.GridAlpha = 0.3;
ax1.FontSize = 22;
% ax1.XLim = [0, 1];
% ax1.YLim = [-50, -30];
% ax1.YLim = [30, 50];
% xscale(ax1, "log");

% Convert negative numbers to LaTeX format with proper minus sign
% yticks(ax1, linspace(-50, -30, 5)) % Adjust tick locations
% yticks(ax1, linspace(30, 50, 5)) % Adjust tick locations

% yticklabels(strrep(string(yticks), '-', char(8722)))
% ytickformat('%.1f')
xtickformat('%.1f')
xtickangle(ax1, 0);
% Enable LaTeX interpreter
% ax1.TickLabelInterpreter = 'latex';
ax1.YAxis.TickLabelInterpreter = 'none';

%% Miscellenous Function
function [new_data] = data_consistency(old_data, current_data)
    if  ~(size(old_data, 1) == size(current_data, 1))
        num_row_diff = size(old_data, 1) - size(current_data, 1);
        if num_row_diff>0
            current_data = vertcat(current_data, NaN(num_row_diff, size(current_data, 2)));
        elseif num_row_diff<0
            old_data = vertcat(old_data, NaN(abs(num_row_diff), size(old_data, 2)));
        end
    end
    new_data = horzcat(old_data, current_data);
end

function [] = consistent_color(ax)
    
    % Grab handles of scatter plots
    scatter_handles = findobj(ax, 'Type', 'Scatter');

    % Grab handles of line plots
    line_handles = findobj(ax, 'Type', 'Line');
    % if isempty(line_handles)
    %     line_handles = 0;
    % end

    % Get color order from the axes
    colors = ax.ColorOrder;

    for i=1:length(scatter_handles)
        % Setting the corresponding scatter and line pair to have the same
        % color index
        if ~isempty(scatter_handles)
            set(scatter_handles(i), 'CData', colors(i, :));
        end

        if ~isempty(line_handles)
            set(line_handles(i), 'Color', colors(i, :));
        end
    end
end