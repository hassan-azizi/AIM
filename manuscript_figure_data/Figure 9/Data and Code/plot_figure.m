% function [] = ss()
clc;
clear;
close all;
%% Load Data
CO2_CALF_Data = readmatrix("./data/CALF_20_CO2_298K.csv");
P = CO2_CALF_Data(:, 1)';
q = CO2_CALF_Data(:, 2)';
N = length(P);
q_mean = mean(q);
%% Isotherm Function
iso_model = 'DS-Langmuir';
iso_struc = isotherm_fit_opt(iso_model, q, P);
iso_fun = iso_struc.fun;

%% Fitted Parameters
% method_array = {'Mean Distribution', 'IsoFit', 'RUPTURA', 'IAST++', 'pyIAST'};
method_array = {'IsoFit'};

method_array  = flip(method_array);
% fit_params = [dist_params(1:end-1, 1)';...
%               2.75681, 8.02640e-04, 3.19182, 4.65637e-06;...
%               3.20806, 4.08438e-06, 2.81696, 7.65120e-04;...
%               2.91204, 6.81999e-04, 3.27074, 3.33057e-06;...
%               3.20808, 4.00000e-06, 2.81698, 7.65000e-04;...
%               ];   % R1: UA, R2: AIM IsoFiT, R3: RUPTURA, R4: IASTpp, R5: pyIAST

% fit_params = [dist_params(1:end-1, 1)';...
%               2.81704, 7.65077e-04, 3.20816, 4.08352e-06;...
%               3.20806, 4.08438e-06, 2.81696, 7.65120e-04;...
%               2.91204, 6.81999e-04, 3.27074, 3.33057e-06;...
%               3.20808, 4.00000e-06, 2.81698, 7.65000e-04;...
%               ];   % R1: UA, R2: AIM IsoFiT, R3: RUPTURA, R4: IASTpp, R5: pyIAST

fit_params = [2.81696, 7.65120e-04, 3.20806, 4.08438e-06];   % ISoFIT

fit_params = flip(fit_params, 1);
M = size(fit_params, 2);    % Number of params in DS-Langmuir Function
%
%% RMSE Calculation
RMSE_values = zeros(length(method_array), 1);
r2_values = zeros(length(method_array), 1);
q_pred = zeros(length(method_array), N);

for i=1:length(method_array)
    y = iso_fun(fit_params(i, :), P);
    SSE = sum((q-y).^2);

    RMSE_values(i) = sqrt(SSE/(N-M));
    r2_values(i) = 1 - SSE/(sum((q-q_mean).^2));
    q_pred(i, :) = y;
end

%% Plots
fig1 = figure(1);
fig1.Position = [100, 50, 850, 750];
% color_array = flip(["red", "blue", "green", "magenta"]);
% color_array = ["red"];
hold on


s_p = scatter(P, q, 'filled');
s_p.SizeData = 110;
s_p.MarkerEdgeColor = 'k';
s_p.MarkerFaceColor = 'k';
s_p.LineWidth = 0.4;
s_p.MarkerFaceAlpha = 0.60; 
% s_p.MarkerFaceColor = color_array(i); 
s_p.LineWidth = 1.0;    

for i=1:length(method_array)
    l_p = plot(P, q_pred(i, :));
    l_p.LineWidth = 2.5;
    % l_p.Color = color_array(i);
end

hold off

ax1 = gca();
xlabel(ax1, "Pressure (Pa)", FontSize=26);

ylabel(ax1, "CO_{2} Uptake (mol/kg)", FontSize=26)
legend(ax1, [{'GCMC data'}, method_array], FontSize=20, Location="best", NumColumns=1);
% legend(ax1, {'IsoFit', 'GCMC'}, FontSize=20, Location="best", NumColumns=1);

grid on
box on
ax1.FontName = 'Deja Vu Sans';
ax1.GridAlpha = 0.3;
ax1.FontSize = 22;
ax1.XLim = [10e-0, 10e5];
ax1.XTick = [10e-0, 10e1, 10e2, 10e3, 10e4, 10e5];
% ax1.YLim = [0, 2.0];
xscale(ax1, "log");
ytickformat('%.1f')
xtickangle(ax1, 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig2 = figure(2);
fig2.Position = [100, 50, 850, 750];
% color_array = ["#1fa720", "#c33919", "m"];
hold on
for i=1:length(method_array)
    s_p = scatter(q, q_pred(i, :), 'filled');
    s_p.SizeData = 110;
    s_p.MarkerEdgeColor = 'k';
    s_p.LineWidth = 0.4;
    s_p.MarkerFaceAlpha = 0.80; 
    % s_p.MarkerFaceColor = color_array(i); 
    s_p.LineWidth = 1.0;   
end

plot(q, q, '-k', 'LineWidth',2.5);

hold off

ax1 = gca();
xlabel(ax1, "Actual Gas Uptake (mol/kg)", FontSize=26);

ylabel(ax1, "Predicted Gas Uptake (mol/kg)", FontSize=26)
legend(ax1, method_array, FontSize=20, Location="best", NumColumns=1);

grid on
box on
ax1.FontName = 'Deja Vu Sans';
ax1.GridAlpha = 0.3;
ax1.FontSize = 22;
% ax1.XLim = [0, 600];
% ax1.YLim = [0, 2.0];
% xscale(ax1, "log");
ytickformat('%.1f')
xtickangle(ax1, 0);

% opt = optimset('TolFun', 1e-10, 'TolX', 1e-10, 'MaxFunEvals', 1e6, 'MaxIter', 1e6);
% s = fminsearch(@res_fun, zeros(1, 4), opt)
% g = 0;
% % Res fun
% function res = res_fun(x)
% if any(x<0)
%     res=1e03;
% else
%     res=sum((iso_fun(x, P)-q).^2);
% end
% end
% end