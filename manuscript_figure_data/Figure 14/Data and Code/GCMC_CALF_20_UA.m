%% Format Screen/Memory
clc;
clear;
%
%% Isotherm Data Processing
% CALF-20 CO2 Data
data_1 = readmatrix("Isotherm Data/CALF-20/Isotherm_273.csv");
data_2 = readmatrix("Isotherm Data/CALF-20/Isotherm_298.csv");
data_3 = readmatrix("Isotherm Data/CALF-20/Isotherm_323.csv");

temperatures = [273, 298, 323];

new_data = data_consistency([], data_1(:, 1:2));
new_data = data_consistency(new_data, data_2(:, 1:2));
new_data = data_consistency(new_data, data_3(:, 1:2));

pressure_data = new_data(:, 1:2:end);
loading_data = new_data(:, 2:2:end);

Pressure_exp = log(rmmissing(reshape(pressure_data, [], 1)));
Loading = rmmissing(reshape(loading_data, [], 1));

len_non_nan_data = sum(~isnan(pressure_data), 1);

temperature_data_m = ones(sum(len_non_nan_data), 1);
k = 0;
for m=1:length(temperatures)    
    temperature_data_m(k+1:k+len_non_nan_data(m), 1) = temperatures(m);
    k = k+len_non_nan_data(m);
end

% Removing zeros loding values and corresponding pressure points
idx = find(Loading==0);
Loading(idx) = [];
Pressure_exp(idx) = [];
temperature_data_m(idx) = [];

Pressure_exp = Pressure_exp(1:2:end);
Loading = Loading(1:2:end);
temperature_data_m = temperature_data_m(1:2:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nominal parameters, obtained from AIM HeatFit
% theta_nom_orig = [-4545.6310, 355.8787, -454.9651, -27.7425, 145.6276, -53.1711, 7.3132, -0.3393,...
%                     21.5519, -0.4417, 0.9066, -0.1107];
% num_a_params = 8;
% param_names = ["a_{0}", "a_{1}", "a_{2}", "a_{3}", "a_{4}", "a_{5}", "a_{6}", "a_{7}",...
%                 "b_{0}", "b_{1}", "b_{2}", "b_{3}"];

% theta_nom_orig = [-4.845673e+03, 9.654406e02, -1.223853e03, 5.838586e02, -1.221194e2, 9.281350,...
%                  22.501361, -1.547888, 1.277603, -0.149109];

% After resdiual function modification
theta_nom_orig = [-4.78E+03, 7.65E+02, -1.11E+03, 5.56E+02, -1.19E+02, 9.079949,...
                    22.294813, -1.031967, 1.055692,-0.124349];
% theta_nom_orig =[-4.69E+03	2.63E+02	-8.55E+02	5.36E+02	-1.23E+02	9.441289...
%                 21.99285	0.688322	0.120564];

theta_nom_orig =[-4.686E+03	2.63E+02	-8.55E+02	5.36E+02	-1.23E+02	9.441289...
                21.99285	0.688322	0.120564];

num_a_params = 6;
param_names = ["a_{0}", "a_{1}", "a_{2}", "a_{3}", "a_{4}", "a_{5}",...
                "b_{0}", "b_{1}", "b_{2}", "b_{3}"];


theta_nominal = abs(theta_nom_orig);

isofun = @(params, loading) Virial_func(params, loading, temperature_data_m, num_a_params);

Pressure_pred = isofun(theta_nom_orig, Loading);


% param_units = ["mol/g", "1/Pa", "mol/g", "1/Pa", "-"];
param_units = repmat("-", length(param_names));

% The standard deviation of error will be treated as a random variable (parameter)
% in analysis and the nominal value is taken from the fitting results using
% nominal parameters.
% error_dist = fitdist((loading_exp-loading_pred)', "Normal");
% error_sigma = error_dist.sigma;
error_nominal = Pressure_exp - Pressure_pred;
error_sigma = std(error_nominal);

% Appending nominal params array
theta_nominal = [theta_nominal, error_sigma];

% Plots
% figure(1);
% hold on
% plot(Loading, Pressure_pred);
% scatter(Loading, Pressure_exp, "filled");
% hold off
% xlabel("Loading (mol/kg)");
% ylabel("ln P (Pa)");
% title("Isotherm Plot using Nominal Parameters");
% 
% figure(2);
% histogram(error_nominal, 50);
% xlabel("Nominal Error (-)");
% ylabel("Count (-)");
% title("Error using Nominal Parameters");
%
%% Bays Analysis and Joint Posterior Distribution
N = 1e08;                           % Number of MC trials
np = length(theta_nominal);         % Number of parameters

% Bound width of MC samples in % of nominal values.
% delta_params = [50, 99.99, 50, 99.99, 30];
% delta_params = [10, 30, 10, 30, 50];
delta_params = zeros(size(theta_nominal)) + 3;
for s=1:np
    if theta_nominal(s) < 10
        delta_params(s) = 10;
    elseif theta_nominal(s) < 100
        delta_params(s) = 3;
    elseif theta_nominal(s) < 1000
        delta_params(s) = 3;
    elseif theta_nominal(s) < 2000
        delta_params(s) = 3;
    else
        delta_params(s) = 3;
    end
end

% delta_params = [10, 10, 10, 10, 10, 10, 10, 10,...
%                 10, 10, 10, 10, 10].* 0.5;

% Create mask for negative model parameters
neg_params = ones(size(theta_nominal));
for i = 1:(np-1)
    if theta_nom_orig(i) < 0
        neg_params(i) = -1;
    end
end

% Defining upper and lower bounds
lb = zeros(1, np);
ub = zeros(1, np);
for i=1:np
    lb(i) = (1-delta_params(i)/100) * theta_nominal(i);
    ub(i) = (1+delta_params(i)/100) * theta_nominal(i);       
end

% Quasi_MC (QMC) sample points generation using Halton Set.
halton_set = haltonset(np);
sobol_set = sobolset(np);
% QMC_samples = net(halton_set, N);
QMC_samples = net(sobol_set, N);
g = getP_D_T(theta_nominal, Pressure_exp, Loading, isofun, neg_params);
% Generating bounded QMC samples
for i=1:N
    for j=1:np
        QMC_samples(i, j) = lb(j) + QMC_samples(i, j)*(ub(j)-lb(j));
    end
end

% Generating MC matrix conatining samples, likelihood, and prior distribution
% for the given samples.
MC_matrix = zeros(N, np+2);
parfor i =1:N
    MC_matrix(i, :) = [QMC_samples(i, :),...
                     getP_D_T(QMC_samples(i, :), Pressure_exp, Loading, isofun, neg_params),...
                     get_P_T(QMC_samples(i, :), theta_nominal)];
end

% Normalizing Constant calculated using QMC integration.
% We are not using area in the integration because the likelihood and prior
%  function returns probabilities which are already normalised to 1. 
PD = (1/size(MC_matrix, 1))*sum(MC_matrix(:, end-1).*MC_matrix(:, end));

% Joint Posterior Distribution
% P_T_D = (MC_matrix(:, end-1) .* MC_matrix(:, end)) ./ PD;
% figure(3);
% histogram(P_T_D, 100);
%
%% Marginal Posterior Distribution
% Number of bins in marginal pos samples
num_bin_marg_pos = 60;
% Matrix containing edge values of bins of samples of marginal posteriors
% for each parameter
marg_pos_bin_edges = zeros(np, num_bin_marg_pos+1); 

for j = 1:np
    % Bins edges are uniformly spaced and are bounded by minimum and 
    % maximum sampled values for each parameter in the QMC sample
	marg_pos_bin_edges(j, :) = linspace(min(QMC_samples(:, j)),...
                                        max(QMC_samples(:, j)),...
                                        size(marg_pos_bin_edges, 2)); 
end

% Sampling marginal posterier distribution for each parameter
marg_pos = zeros(np, num_bin_marg_pos);
integrand_store = zeros(size(marg_pos));
integrand_count = zeros(size(marg_pos));

for s=1:size(MC_matrix, 1)
    for p=1:np
        idx = find(marg_pos_bin_edges(p, :) >= QMC_samples(s, p), 1) - 1;
        if idx > 0
            integrand_store(p, idx) = integrand_store(p, idx)...
                                      + MC_matrix(s, end-1)*MC_matrix(s, end)/PD;
            integrand_count(p, idx) = integrand_count(p, idx) + 1;
        end
    end    
end

marg_pos = integrand_store./integrand_count;
marg_pos(isnan(marg_pos)) = 0;
%
%% Plotting Marginal Priors
% priors_sample = get_priors(theta_nominal, lb, ub, 1000);
% f = figure(4);
% f.Position = [0, 0, 750, 750];
% 
% for p=1:np
%     subplot(3, 2, p);
%     histogram(priors_sample(p, :), num_bin_marg_pos, 'Normalization', 'pdf');
%     ax = gca;
%     xlim(ax, [lb(p), ub(p)]);
%     xlabel(ax, sprintf("%s, (%s)", param_names(p), param_units(p)));
%     ylabel(ax, sprintf("P(%s|D)", param_names(p)));
%     set(ax, 'YTick', []);
% 	% set(ax, 'Fontsize', 15, 'Linewidth', 2);
%     pbaspect([1, 1, 1]);
% end
% sgtitle("Prior Distribution");
%
%% Plotting Marginal Posteriors
f = figure(5);
f.Position = [0, 0, 1000, 1000];

% hist_values = zeros(np, size(marg_pos, 2));     % To be used for fitting distribution
hist_values = marg_pos;
fitted_dist = zeros(np, 3);% Location, scale and shape params 

% Fitted Distributions
for p = 1:np
    % fitted_dist(p, :) = fit_dist(marg_pos_bin_edges(p, :), hist_values(p, :));
end

for p=1:np-1
    subplot(4, 3, p);
    hold on
    % Plot marginal posterior histogram
    h = histogram('BinEdges', marg_pos_bin_edges(p, :),...
              'BinCounts', marg_pos(p, :), 'Normalization', 'pdf');
    % Plot marginal posterior distribution
    % x = linspace(marg_pos_bin_edges(p, 1), marg_pos_bin_edges(p, end), 200);
    % y = get_dist(fitted_dist(p, :), x);
    % plot(x, y, '-r');
    ax = gca;
    low_lim = min([lb(p), ub(p)]);
    high_lim = max([lb(p), ub(p)]);
    xlim(ax, [low_lim, high_lim]);
    xlabel(ax, sprintf("%s, (%s)", param_names(p), param_units(p)));
    ylabel(ax, sprintf("P(%s|D)", param_names(p)));
    % set(ax, 'YTick', []);
    % hold off
	set(ax, 'Fontsize', 12, 'Linewidth', 1.2);
    % pbaspect([1, 1, 1]);
end
subplot(4, 3, 12);
% Q_ads UA
Q_ads_edges = marg_pos_bin_edges(1, :) .* 8.314e-03;
Q_ads_marg_pos = marg_pos(1, :) .* 8.314e-03;
h = histogram('BinEdges', Q_ads_edges,...
              'BinCounts', Q_ads_marg_pos, 'Normalization', 'pdf');
ax = gca;
low_lim = min([lb(1), ub(1)].*8.314e-03);
high_lim = max([lb(1), ub(1)].*8.314e-03);
xlim(ax, [low_lim, high_lim]);
xlabel(ax, sprintf("%s, (%s)", "Q_{ads}", "kJ/mole"));
ylabel(ax, sprintf("P(%s|D)", "Q_{ads}"));
set(ax, 'YTick', []);
set(ax, 'Fontsize', 12, 'Linewidth', 1.2);
% pbaspect([1, 1, 1]);
sgtitle("Posterior Distribution");

%% Miscellenous Function
% Likelihood function P(Data/theta). Return probability of how likely is it
% to match the experimental data given the set of values of parameters.
function P_D_T = getP_D_T(theta_vec, P_exp, Q_exp, iso_fun, neg_params)
    N_p = length(Q_exp);
    theta_vec = theta_vec.*neg_params;
    P_pred = iso_fun(theta_vec(1:end-1), Q_exp);
    sigma = theta_vec(end);
    SSE = sum((P_exp-P_pred).^2, "all");
    multiplier = 1/(sigma*sqrt(2*pi))^N_p;
    % multiplier = 1/(sigma*sqrt(2*pi));
    exp_part = exp(-SSE/2/sigma^2);

    P_D_T = multiplier.*exp_part;
    if isnan(P_D_T)
        s=100;
    end
end

% Prior Distribution P(theta). Return the probability of how likely is the given
% value of parameter. Here we are using Maximum Entropy Principle.
function P_theta = get_P_T(theta_vec, theta_nominal)
    % Returns the Joint Prior probability for given theta
    P_theta_i = 1./theta_nominal .* exp(-theta_vec./theta_nominal);
    P_theta = prod(P_theta_i);
end

% Prior Marginal samples
function Priors = get_priors(theta_nominal, lb, ub, num_samples)
    np = length(theta_nominal);
    
    Priors = zeros(np, num_samples);

    % Rejection sampling for each parameter
    for i = 1:length(theta_nominal)
        mu = theta_nominal(i);
        L = lb(i);
        U = ub(i);
        count = 0;
        while count < num_samples
            x = exprnd(mu);  % Sample from exponential distribution
            if x >= L && x <= U
                count = count + 1;
                Priors(i, count) = x;
            end
        end
    end
end

% Function for skewed normal distribution
function y = get_dist(params, x)
    % Unpack parameter values	
    e = params(1); % Location parameter
    w = params(2); % Scale parameter
    a = params(3); % Shape parameter
    
    % Caluclate vector of skew-normal PDF values (f)
    phi_pdf = (1/sqrt(2*pi))*exp(-((x-e).^2)/2/(w^2));
    phi_cdf = 0.5*(1 + erf(a*(x-e)/w/sqrt(2)));
    y = 2*phi_pdf.*phi_cdf/w;
end