%% Data Loading
addpath("data\")
data_GUI_CO2 = readmatrix("CO2_data.csv");
data_GUI_N2 = readmatrix("N2_data.csv");
data_GUI_He = readmatrix("Helium.csv");
time_GUI = data_GUI_CO2(:, 1);

%% Composition Plots
fig1 = figure(1);
fig1.Position = [100, 50, 850, 750];

hold on
q = plot(time_GUI, data_GUI_N2(:, end), '-b',...
         time_GUI, data_GUI_He(:, end), '-m',...
         time_GUI, data_GUI_CO2(:, end), '-r');

hold off

for i=1:3
    q(i).LineWidth = 2.5;
end

ax1 = gca();

grid on
ax1.FontName = 'Deja Vu Sans';
ax1.FontSize = 22;
ax1.GridAlpha = 0.3;

xlabel(ax1, "Time (s)", "FontSize",26);
ylabel(ax1, "Mole Fraction", "FontSize",26)

ax1.XLim = [0, 5000];

legend(ax1, {'N_2', 'He', 'CO_2'}, NumColumns=1, FontSize=20);


ytickformat('%.1f')
box on
xtickangle(ax1, 0);