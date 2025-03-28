%% GCMC DATA
rmpath("data\");
addpath("data\");
data_Xe = readmatrix("SBMOF_1_xenon_298K.csv");
data_Kr = readmatrix("SBMOF_1_krypton_298K.csv");

GCMC_Xe_loading = data_Xe(:, 2);
GCMC_Kr_loading = data_Kr(:, 2);
pressure = data_Xe(:, 1);

%% Loading calculation by Langmuir model 

Xe = Isotherm("DS-Langmuir", [1.46, 9.68E-04, 0, 1, 0, 0, 0, 1]);
Kr = Isotherm("DS-Langmuir", [1.47, 2.92E-05, 0, 1, 0, 0, 0, 1]);
temperature = 298 .* ones(size(pressure));

loading_CO2 = Xe.pure_loading(pressure, temperature);
loading_CH4 = Kr.pure_loading(pressure, temperature);

%% PLOTS
fig1 = figure(2);
fig1.Position = [100, 50, 850, 750];

% p = semilogx(pressure/1e05, loading_CO2, '-r', pressure/1e05, loading_CH4, '-b');
hold on
p = scatter(pressure/1e05, GCMC_Xe_loading, 'cyan', 'filled');
q = loglog(pressure/1e05, loading_CO2, '-c');


for i=1:1
    q(i).LineWidth = 2.5;
    q(i).Color = "#1fa720";
    p(i).SizeData = 90;
    p(i).MarkerEdgeColor = 'k';
    p(i).MarkerFaceColor = "#1fa720";
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
p = scatter(pressure/1e05, GCMC_Kr_loading, 'filled');
q = loglog(pressure/1e05, loading_CH4, '-m');


% p = loglog(pressure/1e05, loading_CH4, '-b');
for i=1:1
    q(i).LineWidth = 2.5;
    q(i).Color = "#c33919";
    p(i).SizeData = 90;
    p(i).MarkerEdgeColor = 'k';
    p(i).MarkerFaceColor = "#c33919";
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

color_array = ["#1fa720", "#c33919"];
hold on
p = [scatter(pressure/1e05, GCMC_Xe_loading, 'cyan', 'filled'), ...
    scatter(pressure/1e05, GCMC_Kr_loading, 'magenta', 'filled')];

q = loglog(pressure/1e05, loading_CO2, '-c', pressure/1e05, loading_CH4, '-m');

for i=1:2
    q(i).LineWidth = 2.5;
    q(i).Color = color_array(i);
    p(i).SizeData = 90;
    p(i).MarkerEdgeColor = 'k';
    p(i).MarkerFaceColor = color_array(i);
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
legend(ax1, {'GCMC data Xe', 'GCMC data Kr', 'SSL fit', 'SSL fit'}, NumColumns=1, FontSize=20);

ax1.XLim = [10^-05, 10^01];
ax1.XTick = [10^-05 10^-04 10^-03 10^-02 10^-01 10^00 10^01];

ytickformat('%.1f')
box on
xtickangle(ax1, 0);