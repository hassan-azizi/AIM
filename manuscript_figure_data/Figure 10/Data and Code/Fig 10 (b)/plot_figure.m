%% GCMC DATA
rmpath("data\");
addpath("data\");
data_CO2 = readmatrix("CALF_20_CO2_298K.csv");
data_CH4 = readmatrix("CALF_20_methane_298K.csv");

GCMC_CO2_loading = data_CO2(:, 2);
GCMC_CH4_loading = data_CH4(:, 2);
pressure = data_CO2(:, 1);

CO2 = Isotherm("DS-Langmuir", [2.82, 7.65E-04, 0, 1, 3.21, 4.08E-06, 0, 1]);
CH4 = Isotherm("DS-Langmuir", [2.87, 2.08E-05, 0, 1, 0, 0, 0, 1]);
temperature = 298 .* ones(size(pressure));
%% PURE LOADING
loading_CO2 = CO2.pure_loading(pressure, temperature);
loading_CH4 = CH4.pure_loading(pressure, temperature);

%% PLOTS
fig1 = figure(2);
fig1.Position = [100, 50, 850, 750];

% p = semilogx(pressure/1e05, loading_CO2, '-r', pressure/1e05, loading_CH4, '-b');
hold on
p = scatter(pressure/1e05, GCMC_CO2_loading, 'red', 'filled');
q = loglog(pressure/1e05, loading_CO2, '-r');


for i=1:1
    q(i).LineWidth = 2.5;
    p(i).SizeData = 90;
    p(i).MarkerEdgeColor = 'k';
    p(i).LineWidth = 0.6;
    p(i).MarkerFaceAlpha = 0.60; 
end

hold off
ax1 = gca();
set(ax1,'xscale','log');
% set(ax1,'yscale','log');
grid on
ax1.FontName = 'Deja Vu Sans';
ax1.FontSize = 22;
ax1.GridAlpha = 0.3;

xlabel(ax1, "Pressure (bar)", "FontSize",26);
ylabel(ax1, "Uptake (mol/kg)", "FontSize",26);
legend(ax1, {'GCMC Data', 'DS-Langmuir fit'}, NumColumns=1, FontSize=20);

ax1.XLim = [10^-05, 10^01];
ax1.XTick = [10^-05 10^-04 10^-03 10^-02 10^-01 10^00 10^01];

ytickformat('%.1f')
box on
xtickangle(ax1, 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig1 = figure(3);
fig1.Position = [100, 50, 850, 750];

% p = semilogx(pressure/1e05, loading_CO2, '-r', pressure/1e05, loading_CH4, '-b');
hold on
p = scatter(pressure/1e05, GCMC_CH4_loading, 'blue', 'filled');
q = loglog(pressure/1e05, loading_CH4, '-b');


% p = loglog(pressure/1e05, loading_CH4, '-b');
for i=1:1
    q(i).LineWidth = 2.5;
    p(i).SizeData = 90;
    p(i).MarkerEdgeColor = 'k';
    p(i).LineWidth = 0.6;
    p(i).MarkerFaceAlpha = 0.60; 
end

hold off
ax1 = gca();
set(ax1,'xscale','log');
% set(ax1,'yscale','log');

grid on
ax1.FontName = 'Deja Vu Sans';
ax1.FontSize = 22;
ax1.GridAlpha = 0.3;

xlabel(ax1, "Pressure (bar)", "FontSize",26);
ylabel(ax1, "Uptake (mol/kg)", "FontSize",26);
legend(ax1, {'GCMC Data', 'DS-Langmuir fit'}, NumColumns=1, FontSize=20);

ax1.XLim = [10^-05, 10^01];
ax1.XTick = [10^-05 10^-04 10^-03 10^-02 10^-01 10^00 10^01];

ytickformat('%.1f')
box on
xtickangle(ax1, 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig1 = figure(4);
fig1.Position = [100, 50, 850, 750];


hold on
p = [scatter(pressure/1e05, GCMC_CO2_loading, 'red', 'filled'), ...
    scatter(pressure/1e05, GCMC_CH4_loading, 'blue', 'filled')];

q = loglog(pressure/1e05, loading_CO2, '-r', pressure/1e05, loading_CH4, '-b');

for i=1:2
    q(i).LineWidth = 2.5;
    p(i).SizeData = 90;
    p(i).MarkerEdgeColor = 'k';
    p(i).LineWidth = 0.6;
    p(i).MarkerFaceAlpha = 0.60; 
end

hold off
ax1 = gca();
set(ax1,'xscale','log');
% set(ax1,'yscale','log');

grid on
ax1.FontName = 'Deja Vu Sans';
ax1.FontSize = 22;
ax1.GridAlpha = 0.3;

xlabel(ax1, "Pressure (bar)", "FontSize",26);
ylabel(ax1, "Uptake (mol/kg)", "FontSize",26);
legend(ax1, {'GCMC data CO_2', 'GCMC data CH_4', 'DSL fit', 'SSL fit'}, NumColumns=1, FontSize=20);

ax1.XLim = [10^-05, 10^01];
ax1.XTick = [10^-05 10^-04 10^-03 10^-02 10^-01 10^00 10^01];

ytickformat('%.1f')
box on
xtickangle(ax1, 0);