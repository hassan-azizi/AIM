clear;
rmpath("data\");
addpath("data\");

%% BreakLab DATA Ext Lang
data_CO2 = readmatrix("CO2_Lang.csv");
data_CH4 = readmatrix("CH4_Lang.csv");
data_helium = readmatrix("Helium_Lang.csv");
MATLAB_CO2_Lang = data_CO2(:, end)./0.05;
MATLAB_CH4_Lang = data_CH4(:, end)./0.05;
MATLAB_helium_Lang = data_helium(:, end)./0.90;
MATLAB_time_Lang = data_helium(:, 1);

%% BreakLab DATA IAST
data_CO2 = readmatrix("CO2_IAST.csv");
data_CH4 = readmatrix("CH4_IAST.csv");
data_helium = readmatrix("Helium_iast.csv");
MATLAB_CO2_IAST = data_CO2(:, end)./0.05;
MATLAB_CH4_IAST = data_CH4(:, end)./0.05;
MATLAB_helium_IAST = data_helium(:, end)./0.90;
MATLAB_time_IAST = data_helium(:, 1);

%% PLOTS
fig2 = figure(2);
fig2.Position = [100, 50, 850, 750];

p = plot(MATLAB_time_IAST, MATLAB_CO2_IAST, '-r', MATLAB_time_IAST, MATLAB_CH4_IAST, '-b', MATLAB_time_IAST, MATLAB_helium_IAST);
for i=1:3
    p(i).LineWidth = 2.5;
    p(3).Color = "m";
end

hold on

q = [scatter(MATLAB_time_Lang(1:30:length(MATLAB_time_Lang(:, 1))), MATLAB_CO2_Lang(1:30:length(MATLAB_time_Lang(:, 1))), 'red', 'filled'),...
         scatter(MATLAB_time_Lang(1:30:length(MATLAB_time_Lang(:, 1))), MATLAB_CH4_Lang(1:30:length(MATLAB_time_Lang(:, 1))), 'blue', 'filled'),...
         scatter(MATLAB_time_Lang(1:30:length(MATLAB_time_Lang(:, 1))), MATLAB_helium_Lang(1:30:length(MATLAB_time_Lang(:, 1))), 'filled')];
for i=1:3
        q(i).SizeData = 90;
        q(i).MarkerEdgeColor = 'k';
        q(i).LineWidth = 0.6;
        q(i).MarkerFaceAlpha = 0.60; 
        q(i).LineWidth = 1.0;
        q(3).MarkerFaceColor = "m";
end
hold off

ax2 = gca();
xlabel(ax2, "Time (s)", FontSize=26);

ylabel(ax2, "Composition y_{i}/y_{io}", FontSize=26)
legend(ax2, {'CO_2 _{IAST}', 'CH_4 _{IAST}', 'He _{IAST}', 'CO_2 _{EDSL}', 'CH_4 _{EDSL}', 'He _{EDSL}'}, ...
     FontSize=20, Location="best", NumColumns=1);

grid on
box on
ax2.FontName = 'Deja Vu Sans';
ax2.GridAlpha = 0.3;
ax2.FontSize = 22;
% ax2.XLim = [0, 600];
ax2.YLim = [0, 1.8];
ytickformat('%.1f')
xtickangle(ax2, 0);