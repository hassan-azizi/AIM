clear;
rmpath("data\");
addpath("data\");
%% RUPTURA DATA
data_CO2 = readmatrix("CO2_RUPTURA.data", "FileType","text");
data_CH4 = readmatrix("CH4_RUPTURA.data", "FileType","text");
data_helium = readmatrix("Helium_RUPTURA.data", "FileType","text");
RUPTURA_CO2 = data_CO2(:, end);
RUPTURA_CH4 = data_CH4(:, end);
RUPTURA_helium = data_helium(:, end);
RUPTURA_time = data_helium(:, 2) .* 60;

%% BreakLab DATA IAST
data_CO2 = readmatrix("CO2_IAST.csv");
data_CH4 = readmatrix("CH4_IAST.csv");
data_helium = readmatrix("Helium_iast.csv");
MATLAB_CO2_IAST = data_CO2(:, end)./0.05;
MATLAB_CH4_IAST = data_CH4(:, end)./0.05;
MATLAB_helium_IAST = data_helium(:, end)./0.90;
MATLAB_time_IAST = data_helium(:, 1);

%% PLOTS
fig1 = figure(1);
fig1.Position = [100, 50, 850, 750];


p = plot(MATLAB_time_IAST, MATLAB_CO2_IAST, '-r', MATLAB_time_IAST, MATLAB_CH4_IAST, '-b', MATLAB_time_IAST, MATLAB_helium_IAST');
for i=1:3
    p(i).LineWidth = 2.5;
    p(3).Color = "m";
end

hold on
q = [scatter(RUPTURA_time(1:20:length(RUPTURA_time(:, 1))), RUPTURA_CO2(1:20:length(RUPTURA_time(:, 1))), 'red', 'filled'),...
    scatter(RUPTURA_time(1:20:length(RUPTURA_time(:, 1))), RUPTURA_CH4(1:20:length(RUPTURA_time(:, 1))), 'blue', 'filled'),...
    scatter(RUPTURA_time(1:20:length(RUPTURA_time(:, 1))), RUPTURA_helium(1:20:length(RUPTURA_time(:, 1))), 'filled')];
for i=1:3
        q(i).SizeData = 90;
        q(i).MarkerEdgeColor = 'k';
        q(i).LineWidth = 0.6;
        q(i).MarkerFaceAlpha = 0.60; 
        q(i).LineWidth = 1.0;
        q(3).MarkerFaceColor = "m";
        
end

ax1 = gca();
xlabel(ax1, "Time (s)", FontSize=26);

ylabel(ax1, "Composition y_{i}/y_{io}", FontSize=26)
legend(ax1, {'CO_2 _{BreakLab}', 'CH_4 _{BreakLab}', 'He _{BreakLab}', 'CO_2 _{RUPTURA}', 'CH_4 _{RUPTURA}', 'He _{RUPTURA}'}, ...
     FontSize=20, Location="best", NumColumns=1);
% legend(ax1, {'Xe IAST', 'Kr IAST', 'He IAST', 'Xe EDSL', 'Kr EDSL', 'He EDSL'}, ...
%      FontSize=20, Location="best");

grid on
box on
ax1.FontName = 'Deja Vu Sans';
ax1.GridAlpha = 0.3;
ax1.FontSize = 22;
% ax1.XLim = [0, 600];
ax1.YLim = [0, 1.8];
ytickformat('%.1f')
xtickangle(ax1, 0);