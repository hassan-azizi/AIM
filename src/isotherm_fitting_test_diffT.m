% addpath("Isotherm Data/");
% MOF-801 CO2 Data
% data_1 = readmatrix("Isotherm Data\CO2_MOF-801_data_273.csv");
% data_2 = readmatrix("Isotherm Data\CO2_MOF-801_data_298.csv");
% data_3 = readmatrix("Isotherm Data\CO2_MOF-801_data_323.csv");
% temperature_data = [273, 298, 323];

% C-PAF H2 Data
% data_1 = readmatrix("Isotherm Data/H2_CPAF_77K.csv");
% data_2 = readmatrix("Isotherm Data/H2_CPAF_87K.csv");
% data_3 = readmatrix("Isotherm Data/H2_CPAF_160K.csv");
% data_4 = readmatrix("Isotherm Data/H2_CPAF_233K.csv");
% data_5 = readmatrix("Isotherm Data/H2_CPAF_253K.csv");
% data_6 = readmatrix("Isotherm Data/H2_CPAF_273K.csv");
% data_7 = readmatrix("Isotherm Data/H2_CPAF_298K.csv");
% temperature_data = [77, 87, 160, 233, 253, 273, 298];
% temperature_data = [77, 87, 160];

% C-PAF mod Data
data_1 = readmatrix("Isotherm Data/H2_Data_mod/H2_CPAF_77K.csv");
data_2 = readmatrix("Isotherm Data/H2_Data_mod/H2_CPAF_87K.csv");
data_3 = readmatrix("Isotherm Data/H2_Data_mod/H2_CPAF_160K.csv");
data_4 = readmatrix("Isotherm Data/H2_Data_mod/H2_CPAF_233K.csv");
data_5 = readmatrix("Isotherm Data/H2_Data_mod/H2_CPAF_253K.csv");
data_6 = readmatrix("Isotherm Data/H2_Data_mod/H2_CPAF_273K.csv");
data_7 = readmatrix("Isotherm Data/H2_Data_mod/H2_CPAF_298K.csv");
temperature_data = [77, 87, 160, 233, 253, 273, 298];

% Si-PAF H2 Data
% data_1 = readmatrix("Isotherm Data\H2_SiPAF_77K.csv");
% data_2 = readmatrix("Isotherm Data\H2_SiPAF_87K.csv");
% data_3 = readmatrix("Isotherm Data\H2_SiPAF_233K.csv");
% temperature_data = [77, 87, 233];

% MIL-53 CO2 Data
% data_1 = readmatrix("Isotherm Data\CO2_MIL-53_195.csv");
% data_2 = readmatrix("Isotherm Data\CO2_MIL-53_223.csv");
% data_3 = readmatrix("Isotherm Data\CO2_MIL-53_273.csv");
% temperature_data = [195, 223, 273];

new_data = data_cosistency([], data_1);
new_data = data_cosistency(new_data, data_2);
new_data = data_cosistency(new_data, data_3);
new_data = data_cosistency(new_data, data_4);
new_data = data_cosistency(new_data, data_5);
new_data = data_cosistency(new_data, data_6);
new_data = data_cosistency(new_data, data_7);

algorithm_id = 'auto';
pressure_data = new_data(:, 1:2:end)./1e03;
loading_data = new_data(:, 2:2:end);

ID_T_ref = 2;
weight_assignment = 'Uniform';
num_params = 8;
num_a_params = 5;
guess_matrix = zeros(num_params, 1);
lb_matrix = -Inf(num_params, 1);
ub_matrix = Inf(num_params, 1);

% % isotherm_struc = Isotherm_fitting(isotherm_model, pressure_data, loading_data, T, 'Uniform', 'true', algorithm_id)
% [dH, dH_RMSE, theta_iso_array, theta_RMSE_array, Isotherm_struc_ref] = Isotherm_fitting_diffT(isotherm_model, pressure_data, loading_data,...
%                                                                            temperature_data, ID_T_ref, algorithm_id, ...
%                                                                            weight_assignment, 'True', 'True');

[dH, dH_inf, Isotherm_struc_ref] = Virial_fitting(pressure_data, loading_data, temperature_data, guess_matrix, lb_matrix, ub_matrix, num_a_params,...
                                                                                weight_assignment, 'True', 1);

function [new_data] = data_cosistency(old_data, current_data)
    if  ~(size(old_data, 1) == size(current_data, 1))
        num_row_diff = size(old_data, 1) - size(current_data, 1);
        if num_row_diff>0
            current_data = vertcat(current_data, NaN(num_row_diff, size(current_data, 2)));
        elseif num_row_diff<0
            old_data = vertcat(old_data, NaN(abs(num_row_diff), size(old_data, 2)));
        end
    end
    new_data = horzcat(old_data, current_data);
end