%% Data Loading
addpath("data\")

data_GUI_T = readmatrix("Temperature.csv");
time_GUI = data_GUI_T(:, 1);
data_GUI_T = data_GUI_T(:, end);
data_GUI_CO2 = readmatrix("CO2_data.csv");
% time_GUI = data_GUI_CO2(:, 1);

%% Temperature Plots
fig1 = figure(2);
fig1.Position = [100, 50, 850, 750];

hold on

q = plot(time_GUI, data_GUI_T(:, end), "Color","#D95319");
hold off

q.LineWidth = 2.5;

ax1 = gca();

grid on
ax1.FontName = 'Deja Vu Sans';
ax1.FontSize = 22;
ax1.GridAlpha = 0.3;

xlabel(ax1, "Time (s)", "FontSize",26);
ylabel(ax1, "Temperature (^oC)", "FontSize",26);

ax1.XLim = [0, 5000];


ytickformat('%.0f')
box on
xtickangle(ax1, 0);