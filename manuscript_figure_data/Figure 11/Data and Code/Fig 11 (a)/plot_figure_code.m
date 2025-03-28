%% RUPTURA DATA
clear;
rmpath("data\");
addpath(".\data\")
data_Xe = readmatrix("Xe_RUPTURA.data", "FileType","text");
data_Kr = readmatrix("Kr_RUPTURA.data", "FileType","text");
data_helium = readmatrix("Helium_RUPTURA.data", "FileType","text");
RUPTURA_Xe = data_Xe(:, end);
RUPTURA_Kr = data_Kr(:, end);
RUPTURA_helium = data_helium(:, end);
RUPTURA_time = data_helium(:, 2) .* 60;

%% BreakLab DATA IAST
data_Xe = readmatrix("Xe_IAST.csv");
data_Kr = readmatrix("Kr_IAST.csv");
data_helium = readmatrix("He_IAST.csv");
MATLAB_Xe_IAST = data_Xe(:, end)./0.05;
MATLAB_Kr_IAST = data_Kr(:, end)./0.05;
MATLAB_helium_IAST = data_helium(:, end)./0.90;
MATLAB_time_IAST = data_helium(:, 1);

%% PLOTS
fig1 = figure(1);
fig1.Position = [100, 50, 850, 750];

color_array = ["#1fa720", "#c33919", "m"];
p = plot(MATLAB_time_IAST, MATLAB_Xe_IAST, '-c', MATLAB_time_IAST, MATLAB_Kr_IAST, '-m', MATLAB_time_IAST, MATLAB_helium_IAST);
for i=1:3
    p(i).LineWidth = 2.5;
    p(i).Color = color_array(i);
end

hold on
q = [scatter(RUPTURA_time(1:10:length(RUPTURA_time(:, 1))), RUPTURA_Xe(1:10:length(RUPTURA_time(:, 1))), 'filled'),...
    scatter(RUPTURA_time(1:10:length(RUPTURA_time(:, 1))), RUPTURA_Kr(1:10:length(RUPTURA_time(:, 1))), 'filled'),...
    scatter(RUPTURA_time(1:10:length(RUPTURA_time(:, 1))), RUPTURA_helium(1:10:length(RUPTURA_time(:, 1))), 'filled')];
for i=1:3
        q(i).SizeData = 90;
        q(i).MarkerEdgeColor = 'k';
        q(i).LineWidth = 0.6;
        q(i).MarkerFaceAlpha = 0.60; 
        q(i).MarkerFaceColor = color_array(i); 
        q(i).LineWidth = 1.0;   
end

ax1 = gca();
xlabel(ax1, "Time (s)", FontSize=26);

ylabel(ax1, "Composition y_{i}/y_{io}", FontSize=26)
legend(ax1, {'Xe_{BreakLab}', 'Kr_{BreakLab}', 'He_{BreakLab}', 'Xe_{RUPTURA}', 'Kr_{RUPTURA}', 'He_{RUPTURA}'}, ...
     FontSize=20, Location="best", NumColumns=1);
% legend(ax1, {'Xe IAST', 'Kr IAST', 'He IAST', 'Xe EDSL', 'Kr EDSL', 'He EDSL'}, ...
%      FontSize=20, Location="best");

grid on
box on
ax1.FontName = 'Deja Vu Sans';
ax1.GridAlpha = 0.3;
ax1.FontSize = 22;
ax1.XLim = [0, 600];
ax1.YLim = [0, 2.0];
ytickformat('%.1f')
xtickangle(ax1, 0);