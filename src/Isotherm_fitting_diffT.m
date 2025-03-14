function [dH, dH_RMSE, theta_iso_array, theta_RMSE_array, Isotherm_struc_ref_2] = Isotherm_fitting_diffT(isotherm_model, ...
                            pressure_data, loading_data, temperature_data, ID_T_ref, guess_matrix, lb_matrix, ub_matrix, ...
                            algorithm_id, weight_assignment, plot_diffT_curves, plot_theta_cuvre, iso_models_available, use_multi)
    %% Function Variables
    % theta_ub = 1000;
    % theta_lb = 0;
    % theta_guess = 1;

    % dH_ub = 50e3;
    % dH_lb = -50e3;
    % dH_guess = 0;

    dH_guess = guess_matrix(end-1, 1);
    dH_lb = lb_matrix(end-1, 1);
    dH_ub = ub_matrix(end-1, 1);

    theta_guess = guess_matrix(end, 1);
    theta_lb = lb_matrix(end, 1);
    theta_ub = ub_matrix(end, 1);
    
    % Reference variables
    pressure_ref = rmmissing(pressure_data(:, ID_T_ref));
    loading_ref = rmmissing(loading_data(:, ID_T_ref));
    T_ref = temperature_data(ID_T_ref);
    
    % Checking that the isotherm data for all the temperature has been
    % provided.
    if ~(size(loading_data, 2) == length(temperature_data))
        error('Isotherm data corresponding to all the temperature values is required...')
    end

    % Checking if all fitting algorithms are to be used
    algorithms_available = {'trust-region-reflective', 'levenberg-marquardt', 'interior-point'};
    if strcmpi(algorithm_id, 'auto')
        fitting_algo_mode = 1;
        fitting_algorithm_array = algorithms_available;
    else
        fitting_algo_mode = 0;
        fitting_algorithm_array = algorithms_available(str2double(algorithm_id)); 
    end
    
    % Handling multistart flag if not specified in input arguments
    if nargin < 14
        use_multi = 0;    
        % use_multi = 1;
    end
    %% Isotherm function fitted for reference temperature
    Isotherm_struc_ref = Isotherm_fitting(isotherm_model, pressure_ref, loading_ref, guess_matrix, lb_matrix, ub_matrix, weight_assignment, 'false', 'auto', iso_models_available, 1.0, use_multi);
    
    %% Predict theta
    theta_iso_array = ones(1, length(temperature_data));
    theta_RMSE_array = zeros(1, length(temperature_data));

    fitting_quantity = 'theta';
    T_array = temperature_data;
    for i = 1:length(temperature_data)
        current_T = temperature_data(i);
        if current_T == T_ref
            continue
        end
        current_pressure = rmmissing(pressure_data(:, i));
        current_loading = rmmissing(loading_data(:, i));
        [final_fit, final_RMSE, ~, ~] = fit_diff_algo(fitting_algo_mode, fitting_algorithm_array, fitting_quantity);
        theta_iso_array(i) = final_fit;
        theta_RMSE_array(i) = final_RMSE;
    end

    %% Predict dH
    fitting_quantity = 'dH';
    % T_array = temperature_data;
    [final_fit, final_RMSE, ~, ~] = fit_diff_algo(fitting_algo_mode, fitting_algorithm_array, fitting_quantity);
    dH = final_fit;
    dH_RMSE = final_RMSE;
    %
    
    %% Fitting the isotherm parameters again using reference temperature function using normalized pressure and loadings
    norm_pressure = reshape(pressure_data.*theta_iso_array, [], 1);
    norm_loadings = reshape(loading_data, [], 1);    
    
    norm_pressure = rmmissing(norm_pressure);
    norm_loadings = rmmissing(norm_loadings);

    % Using the best isotherm model from reference temperature fitting results
    isotherm_model_ref = Isotherm_struc_ref.name;
    Isotherm_struc_ref_2 = Isotherm_fitting(isotherm_model_ref, norm_pressure, norm_loadings, guess_matrix, lb_matrix, ub_matrix, weight_assignment, 'false', 'auto', iso_models_available, 1.0, use_multi);
    
    % Appending dH to isotherm structure array
    Isotherm_struc_ref_2.dH = dH;
    Isotherm_struc_ref_2.dH_RMSE = dH_RMSE;
    Isotherm_struc_ref_2.T_ref = T_ref;
    
    % Get Covarience and dH standard error
    [~, std_error] = get_cov(dH);
    std_error = full(std_error);
    Isotherm_struc_ref_2.dH_std_error = std_error;

    % Appending the RMSE and r2 values for all isotherm model only when user selects auto isotherm model 
    if strcmpi(isotherm_model, 'Auto')
        Isotherm_struc_ref_2.RMSE_all = Isotherm_struc_ref.RMSE_all;
        Isotherm_struc_ref_2.r2_all = Isotherm_struc_ref.r2_all;
    end
    %% Plotting
    if plot_diffT_curves
        % Extracting data from isotherm options structure
        isotherm_fun_ref = Isotherm_struc_ref_2.fun;
        params_ref = Isotherm_struc_ref_2.fitted_params;
        predicted_loadings = zeros(size(pressure_data, 1), size(temperature_data, 2));
        % isotherm_fun_ref = Isotherm_struc_ref.fun;
        % params_ref = Isotherm_struc_ref.fitted_params;
    
        colors = ["blue", "green", "red"];
    
        norm_pressure_data = pressure_data .* exp(dH/8.3144 .* (1./T_array - 1/T_ref));
        
        for k=1:length(T_array)      
            predicted_loadings(:, k) = isotherm_fun_ref(params_ref, norm_pressure_data(:, k)); 
            figure(1)
            scatter(pressure_data(:, k), loading_data(:, k), colors(k));
            hold on
            plot(pressure_data(:, k), predicted_loadings(:, k), colors(k)); 
        end
        hold off
    end

    if plot_theta_cuvre
        figure(2)
        scatter(T_array, theta_iso_array);
        hold on
        plot(T_array, exp(dH/8.3144 .* (1./T_array - 1/T_ref))); 
        hold off
    end
    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                           %% MISCELLENEOUS FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% dH residual Function
    function residual = dH_residual_func(dH, T_ref, T_array, theta_iso_array)        
        % Calculating theta_cc from dH
        theta_cc = exp(-1.0 * dH/8.3144 * (1./T_array - 1/T_ref));       
        % residual = (theta_iso_array - theta_cc).^2;
        residual = (theta_iso_array - theta_cc);
    end
    
    %% Theta residual Function
    function residual = theta_residual_func(theta, isotherm_fun_ref, params_ref, pressure, loading)
        
        % Normalizing the pressure data
        pressure_norm = pressure .* theta;
        loading_pred = isotherm_fun_ref(params_ref, pressure_norm);
        % residual = (loading - loading_pred).^2;
        residual = (loading - loading_pred);
    end
    %% Theta fitting function
    function [fitted_theta_current_algo, resnorm, exitflag, RMSE] = fit_theta(current_options, current_algo)
        
        % Extracting data from isotherm options structure
        isotherm_fun_ref = Isotherm_struc_ref.fun;
        params_ref = Isotherm_struc_ref.fitted_params;

        % Residual function handle
        fit_function = @(x)theta_residual_func(x, isotherm_fun_ref, params_ref, current_pressure, current_loading);
        try
            [fitted_theta_current_algo, resnorm, residual, exitflag] = lsqnonlin(fit_function, theta_guess, theta_lb, theta_ub, current_options);
        
            if exitflag == 0
                error("Maximum Iterations for algorithm %s reached! Terminating the solution...", ...
                        current_algo);
            elseif exitflag < 0
                error("Solver failed to find the isotherm parameters for algorithm %s", ...
                        current_algo);
            end

            % Root Mean Squared Error
            % RMSE = sqrt(mean(residual));   % Think about this definition more, in the denominator there should be a number      
            residual_value = resnorm; 
            RMSE = sqrt(residual_value/length(T_array));

        catch ME1
            disp(ME1.message);
            fitted_theta_current_algo = NaN;
            RMSE = 1e05;
            resnorm = 1e05;
            exitflag = -2;
            % error('%s:\n%s', ME1.identifier, ME1.message);
        end
    end

    %% dH fitting function
    function [fitted_dH_current_algo, resnorm, exitflag, RMSE] = fit_dH(current_options, current_algo)
        
        % Residual function handle
        fit_function = @(x)dH_residual_func(x, T_ref, T_array, theta_iso_array);
        try
            [fitted_dH_current_algo, resnorm, residual, exitflag] = lsqnonlin(fit_function, dH_guess, dH_lb, dH_ub, current_options);
        
            if exitflag == 0
                error("Maximum Iterations for algorithm %s reached! Terminating the solution...", ...
                        current_algo);
            elseif exitflag < 0
                error("Solver failed to find the isotherm parameters for algorithm %s",...
                        current_algo);
            end

            % Root Mean Squared Error
            % residual_value = sum(residual); 
            residual_value = resnorm; 
            RMSE = sqrt(residual_value/length(T_array));   % Think about this definition more, in the denominator there should be a number      
            

        catch ME1
            disp(ME1.message);
            fitted_dH_current_algo = NaN;
            RMSE = 1e05;
            resnorm = 1e05;
            exitflag = -2;
            % error('%s:\n%s', ME1.identifier, ME1.message);
        end
    end
    
    %% Function for fitting theta/dH for different fitting algorithm
    function [final_fit, min_RMSE, exitflag_best, best_algorithm] = fit_diff_algo(fitting_algo_mode, fitting_algorithm_array, fitting_quantity)
        % Fitting using different algorithms
        if fitting_algo_mode  
            RMSE_diff_algo = zeros(1, length(fitting_algorithm_array));
            exit_flag_diff_algo = ones(1, length(fitting_algorithm_array));
            fitted_quantity_diff_algo = zeros(1, length(fitting_algorithm_array));

            for j=1:length(fitting_algorithm_array)
                current_algo = char(fitting_algorithm_array(j));
                % Fitting options
                options = optimset('Display', 'off', 'Diagnostics', 'off', 'MaxIter', 100000,...
                                    'Algorithm', current_algo, 'TolFun', 1e-12, 'TolX', 1e-12,...
                                    'MaxFunEvals', 10000);
                if strcmpi(fitting_quantity, 'theta')
                    [fitted_quantity_current_algo, ~, exitflag, RMSE] = fit_theta(options, current_algo);
                else
                    [fitted_quantity_current_algo, ~, exitflag, RMSE] = fit_dH(options, current_algo);
                end

                RMSE_diff_algo(j) = RMSE;
                exit_flag_diff_algo(j)= exitflag;
                fitted_quantity_diff_algo(j) = fitted_quantity_current_algo;
            end
            
            % Choosing the best algorithm
            [min_RMSE, idxmin] = min(RMSE_diff_algo);
            algo_min_RMSE = char(fitting_algorithm_array(idxmin));
            final_fit = fitted_quantity_diff_algo(idxmin);
            exitflag_best = exit_flag_diff_algo(idxmin);
            best_algorithm = algo_min_RMSE;

        else 
            current_algo = char(fitting_algorithm_array);
            % Fitting options
            options = optimset('Display', 'off', 'Diagnostics', 'off', 'MaxIter', 100000,...
                               'Algorithm', current_algo, 'TolFun', 1e-12, 'TolX', 1e-12,...
                               'MaxFunEvals', 10000);
            if strcmpi(fitting_quantity, 'theta')
                [fitted_quantity_current_algo, ~, exitflag, RMSE] = fit_theta(options, current_algo);
            else
                [fitted_quantity_current_algo, ~, exitflag, RMSE] = fit_dH(options, current_algo);
            end
            
            min_RMSE = RMSE;
            final_fit = fitted_quantity_current_algo;
            exitflag_best = exitflag;
            best_algorithm = current_algo;
        end
    end
    %% Function for estimating covarience matrix and standard error
    function [cov_mat, std_error] = get_cov(fitted_parameters)
        
        % Residual function handle
        fit_function = @(x)dH_residual_func(x, T_ref, T_array, theta_iso_array);
        
        options = optimset('Display', 'off', 'Diagnostics', 'off', 'MaxIter', 100000,...
                                'TolFun', 1e-12, 'TolX', 1e-12,...
                                'MaxFunEvals', 10000);
        try
            % Just a single run of lsqnonlin to get Jacobian
            [~,resnorm, ~,~,~,~,J] = lsqnonlin(fit_function, fitted_parameters, [], [], options);
        
            % Mean Squared Error
            MSE = resnorm/(length(T_array));
            
            % Covarience matrix
            cov_mat = MSE .* inv(J'*J);
            std_error = sqrt(diag(cov_mat));
        catch
            cov_mat = NaN * ones(1, 1);
            std_error = NaN * ones(1, 1);
        end
    end
end