clear;
rmpath("data\");
addpath("data\");
%% BreakLab DATA Ext Lang
data_Xe = readmatrix("Xe_Lang.csv");
data_Kr = readmatrix("Kr_Lang.csv");
data_helium = readmatrix("He_Lang.csv");
MATLAB_CO2_Lang = data_Xe(:, end)./0.05;
MATLAB_CH4_Lang = data_Kr(:, end)./0.05;
MATLAB_helium_Lang = data_helium(:, end)./0.90;
MATLAB_time_Lang = data_helium(:, 1);

%% BreakLab DATA IAST
data_Xe = readmatrix("Xe_IAST.csv");
data_Kr = readmatrix("Kr_IAST.csv");
data_helium = readmatrix("He_IAST.csv");
MATLAB_Xe_IAST = data_Xe(:, end)./0.05;
MATLAB_Kr_IAST = data_Kr(:, end)./0.05;
MATLAB_helium_IAST = data_helium(:, end)./0.90;
MATLAB_time_IAST = data_helium(:, 1);

%% PLOTS
fig2 = figure(2);
fig2.Position = [100, 50, 850, 750];
color_array = ["#1fa720", "#c33919", "m"];

p = plot(MATLAB_time_IAST, MATLAB_Xe_IAST, '-c', MATLAB_time_IAST, MATLAB_Kr_IAST, '-m', MATLAB_time_IAST, MATLAB_helium_IAST);
for i=1:3
    p(i).LineWidth = 2.5;
    p(i).Color = color_array(i);
end
hold on

q = [scatter(MATLAB_time_Lang(1:20:length(MATLAB_time_Lang(:, 1))), MATLAB_CO2_Lang(1:20:length(MATLAB_time_Lang(:, 1))), 'cyan', 'filled'),...
         scatter(MATLAB_time_Lang(1:20:length(MATLAB_time_Lang(:, 1))), MATLAB_CH4_Lang(1:20:length(MATLAB_time_Lang(:, 1))), 'magenta', 'filled'),...
         scatter(MATLAB_time_Lang(1:20:length(MATLAB_time_Lang(:, 1))), MATLAB_helium_Lang(1:20:length(MATLAB_time_Lang(:, 1))), 'filled')];
for i=1:3
        q(i).SizeData = 90;
        q(i).MarkerEdgeColor = 'k';
        q(i).LineWidth = 0.6;
        q(i).MarkerFaceAlpha = 0.60; 
        q(i).LineWidth = 1.0;
        q(i).MarkerFaceColor = color_array(i);
end
hold off

ax2 = gca();
xlabel(ax2, "Time (s)", FontSize=26);

ylabel(ax2, "Composition y_{i}/y_{io}", FontSize=26)
legend(ax2, {'Xe_{IAST}', 'Kr_{IAST}', 'He_{IAST}', 'Xe_{EDSL}', 'Kr_{EDSL}', 'He_{EDSL}'}, ...
     FontSize=20, Location="best", NumColumns=1);

grid on
box on
ax2.FontName = 'Deja Vu Sans';
ax2.GridAlpha = 0.3;
ax2.FontSize = 22;
ax2.XLim = [0, 600];
ax2.YLim = [0, 2.0];
ytickformat('%.1f')
xtickangle(ax2, 0);