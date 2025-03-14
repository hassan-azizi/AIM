function [dH, dH_inf, isotherm_struc] = Virial_fitting(pressure_data,...
                                                loading_data, temperature_data, guess_matrix, lb_matrix, ub_matrix, ...
                                                num_a_params, weight_assignment, plot_res, plot_heat, use_multi)
    %% Data Preparation
    pressure_data_m = rmmissing(reshape(pressure_data, [], 1));
    loading_data_m = rmmissing(reshape(loading_data, [], 1));
    len_non_nan_data = sum(~isnan(pressure_data), 1);
    
    temperature_data_m = ones(sum(len_non_nan_data), 1);
    k = 0;
    for m=1:length(temperature_data)    
        temperature_data_m(k+1:k+len_non_nan_data(m), 1) = temperature_data(m);
        k = k+len_non_nan_data(m);
    end
    
    % Pressure column containe log of pressure values because function
    % returns log of pressure values
    isotherm_data = horzcat(log(pressure_data_m), loading_data_m, temperature_data_m);
    % isotherm_data = horzcat(pressure_data_m, loading_data_m, temperature_data_m);

    % Retaining data which corresponds to only non-zero loading values to
    % avoid error in Virial Fitting.    
    [idx_ret, ~] = find(isotherm_data(:, 2)~=0.0);
    isotherm_data = isotherm_data(idx_ret, :);
    N = size(isotherm_data, 1);
    R = 8.314;  %J/mol

    % Weights assignments correponding to the data points
    switch(weight_assignment)
        case 'Uniform'
            weight_vector = ones(N, 1);
        case 'Biased'
            % Higher weight values for smaller loadings
            % weight_vector = ones(N, 1) .* mean(isotherm_data(: ,2)) ./ isotherm_data(: ,2);
            weight_vector = ones(N, 1) .* 1./(1+isotherm_data(: ,2));
        otherwise
            error("Please use the valid weight assignment method...")
    end

    % Checking if all fitting algorithms are to be used
    algorithms_available = {'trust-region-reflective', 'levenberg-marquardt', 'interior-point'};
    fitting_algo_mode = 1;
    fitting_algorithm_array = algorithms_available;
    
    % Calling Isotherm options
    isotherm_model = 'Virial';
    isotherm_opt = isotherm_fit_opt(isotherm_model, 1);
    
    % Handling multistart flag if not specified in input arguments
    if nargin < 11
        % use_multi = 0;    
        use_multi = 1;    
    end

    % Changing the default values of initial guesses/upper and
    % lower bounds to user specified values. The initial guess will be override if 
    % the user chooses multistart 
    isotherm_opt.num_params = size(guess_matrix, 1);
    isotherm_opt.num_a = num_a_params;
    isotherm_opt.guess = guess_matrix(1:isotherm_opt.num_params, 1)';
    isotherm_opt.lb = lb_matrix(1:isotherm_opt.num_params, 1)';
    isotherm_opt.ub = ub_matrix(1:isotherm_opt.num_params, 1)';
    
    % Initialize and problem and multistart object structures if multistart is used
    if use_multi
        num_trials = 1000;
        multi_problem = createOptimProblem('lsqnonlin');
        multi_object = MultiStart("UseParallel", false, "Display","final",...
                                   "FunctionTolerance",1e-10, "MaxTime", Inf,...
                                   "StartPointsToRun", "all", "XTolerance", 1e-10);
    end

    %% Fitting using lsqnonlin different algorithms and GA
    % [fitted_params_lsq, RMSE_lsq, r2_lsq, exitflag_lsq, ~] = fit_isotherm_diff_algo(isotherm_opt, fitting_algo_mode, fitting_algorithm_array);

    % Using GA to improve the solution
    % [fitted_params_ga, RMSE_ga, r2_ga, exitflag_ga] = run_ga(isotherm_opt, []);
    %
    % isotherm_opt.guess = fitted_params_ga;
    [fitted_params_lsq, RMSE_lsq, r2_lsq, exitflag_lsq, ~] = fit_isotherm_diff_algo(isotherm_opt, fitting_algo_mode, fitting_algorithm_array);

    % if RMSE_ga < RMSE_lsq
    %     best_exitflag = exitflag_ga;
    %     best_r2 = r2_ga;
    %     best_fitted_params = fitted_params_ga;
    %     best_RMSE = RMSE_ga;
    % 
    % else
    %     best_exitflag = exitflag_lsq;
    %     best_r2 = r2_lsq;
    %     best_fitted_params = fitted_params_lsq;
    %     best_RMSE = RMSE_lsq;
    % end
    best_exitflag = exitflag_lsq;
    best_r2 = r2_lsq;
    best_fitted_params = fitted_params_lsq;
    best_RMSE = RMSE_lsq;
    %
    %% Organizing output

    % Get Covarience and standard error
    [~, std_error] = get_cov(best_fitted_params, isotherm_opt);
    std_error = full(std_error);
    
    isotherm_struc.name = isotherm_opt.name;
    isotherm_struc.fun = isotherm_opt.fun;
    isotherm_struc.fitted_params = best_fitted_params;
    isotherm_struc.RMSE = best_RMSE;
    isotherm_struc.r2 = best_r2;
    isotherm_struc.EF = best_exitflag;
    isotherm_struc.T_flag = isotherm_opt.T_flag;
    isotherm_struc.p_sat = isotherm_opt.p_sat;
    isotherm_struc.std_error = std_error;

    if best_r2 <= 0 &&  best_exitflag < 0
        error(['Failed to fit the isotherm. Try another isotherm model or ' ...
                'use better initial guesses...!']);
    end
    
    %% Calculation of Isostearic heat of adsorption
    dH = -1.0 * (-R .* polyval(flip(best_fitted_params(1:num_a_params)), isotherm_data(:, 2)));
    % dH at infinite dilution
    dH_inf = -1.0 * (-R*best_fitted_params(1));
    dH_inf_std_error = R * std_error(1); 
    %
    isotherm_struc.dH = dH;
    isotherm_struc.dH_inf = dH_inf;
    isotherm_struc.dH_inf_std_error = dH_inf_std_error;
    %% Plotting
    if strcmpi(plot_res, "True")
        figure(1);
        cla;
        hold on
        legend_str = {};
        for k=1:length(temperature_data)
            current_loading = rmmissing(loading_data(:, k));
            current_temperature = temperature_data(k);
            pressure_pred = exp(isotherm_struc.fun(isotherm_struc.fitted_params, current_loading, current_temperature, num_a_params));
            
            plot(pressure_pred, current_loading); % Plots for predicted pressures
            scatter(rmmissing(pressure_data(:, k)), current_loading, 'red', 'filled');      % Plots using actual data
            legend_str = [legend_str, {sprintf('Data T = %.2f', temperature_data(k))},...
                            {sprintf('Prediction T = %.2f', temperature_data(k))}];
        end

        hold off
        legend(legend_str);
        xlabel('Pressure (Pa)');
        ylabel('Adsorption Loading (mol/kg)');
    end
    
    if plot_heat
        figure(2);
        cla;
        scatter(isotherm_data(:, 2), dH);
        xlim([0, 5]);
        ylim([3.5, 6.5]);
        xlabel('Adsorption Loading (mol/kg)');
        ylabel('Heat of Adsorption (kJ/mol)');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                %% MISCELLENEOUS FUNCTIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Residual Function
    function residual = residual_func(params, isotherm_fun, pressure, loading, temperatures, num_a_params, weights)    
        pred = isotherm_fun(params, loading, temperatures, num_a_params);
        % pred = exp(isotherm_fun(params, loading, temperatures, num_a_params));
        % residual = weights .* (pressure - pred).^2;
        residual = weights .* (pressure - pred);
    end

    %% Isotherm Fitting Function
    function [fitted_params, resnorm, exitflag, RMSE, r2] = fit_isotherm(isotherm_opt, current_options, current_algo)
        
        % Extracting data from isotherm options structure
        current_model = isotherm_opt.name;
        isotherm_fun = isotherm_opt.fun;
        params_lb = isotherm_opt.lb;
        params_ub = isotherm_opt.ub;
        initial_guess = isotherm_opt.guess;
        num_a = isotherm_opt.num_a;

        % Residual function handle
        fit_function = @(x)residual_func(x, isotherm_fun, isotherm_data(:, 1), ...
                                        isotherm_data(:, 2), isotherm_data(:, 3), num_a, weight_vector);
        try
            if ~use_multi
                [fitted_params, resnorm, ~, exitflag] = lsqnonlin(fit_function, initial_guess, params_lb, params_ub, current_options);
            
            else
                % Update problem structure with new options, if multistart is used
                multi_problem.objective = fit_function;
                multi_problem.x0 = initial_guess;
                multi_problem.lb = params_lb;
                multi_problem.ub = params_ub;
                multi_problem.options = current_options;

                [fitted_params, resnorm, exitflag, output] = run(multi_object, multi_problem, num_trials);
            end

            if exitflag == 0
                error("Maximum Iterations for isotherm model %s and algorithm %s reached! Terminating the solution...", ...
                        current_model, current_algo);
            elseif exitflag < 0
                error("Solver failed to find the isotherm parameters for isotherm model %s and algorithm %s",...
                        current_model, current_algo);
            end

            % Root Mean Square Error
            % residual_value = sum(fit_function(fitted_params));
            % residual_value = sum(residual);
            residual_value = resnorm;
            RMSE = sqrt(residual_value./(N-isotherm_opt.num_params));

            % Coefficient of Determination
            ss_residual = resnorm;
            % ss_residual = sum(residual);
            % ss_residual = sum(fit_function(fitted_params).^2);
            sq_var = sum((isotherm_data(:, 1)-mean(isotherm_data(:, 1), 1)).^2);
            r2 = 1 - ss_residual/sq_var;
            
        catch ME1
            disp(ME1.message);
            fitted_params = NaN(size(initial_guess));
            RMSE = 1e05;
            resnorm = 1e05;
            r2 = 0;
            exitflag = -2;
            % error('%s:\n%s', ME1.identifier, ME1.message);
        end
    end
    
    %% Function for fitting isotherms using different fitting algorithm
    function [fitted_params_current_iso, min_RMSE_current_iso, r2_value_current_iso, exitflag_current_iso, best_algorithm] = fit_isotherm_diff_algo(isotherm_opt, fitting_algo_mode, fitting_algorithm_array)
        % Fitting using different algorithms
        if fitting_algo_mode   
            r2_value_diff_algo = zeros(1, length(fitting_algorithm_array));
            RMSE_diff_algo = zeros(1, length(fitting_algorithm_array));
            exit_flag_diff_algo = ones(1, length(fitting_algorithm_array));
            fitted_param_struc_diff_algo = struct;
            
            for j=1:length(fitting_algorithm_array)
                current_algo = char(fitting_algorithm_array(j));
                % Fitting options
                options = optimset('Display', 'off', 'Diagnostics', 'off', 'MaxIter', 10000000,...
                                   'Algorithm', current_algo, 'TolFun', 1e-12, 'TolX', 1e-12,...
                                   'MaxFunEvals', 10000);
                [fitted_params, ~, exitflag, RMSE, r2_value] = fit_isotherm(isotherm_opt, options, current_algo);
                
                RMSE_diff_algo(j) = RMSE;
                r2_value_diff_algo(j) = r2_value;
                exit_flag_diff_algo(j)= exitflag;
                fitted_param_struc_diff_algo.(strrep(current_algo, '-','_')) = fitted_params;
            end
        
            % Choosing the best algorithm
            [min_RMSE, idxmin] = min(RMSE_diff_algo);
            algo_min_RMSE = char(fitting_algorithm_array(idxmin));
            fitted_params_current_iso = fitted_param_struc_diff_algo.(strrep(algo_min_RMSE, '-','_'));
            min_RMSE_current_iso = min_RMSE;
            exitflag_current_iso = exit_flag_diff_algo(idxmin);
            r2_value_current_iso = r2_value_diff_algo(idxmin);
            best_algorithm = algo_min_RMSE;            
        end
    end
    %
    %% GA optimization function
    %% Residual Function for GA
    function residual = residual_func_ga(params, isotherm_fun, pressure, loading, temperatures, num_a_params, weights)    
        pred = isotherm_fun(params, loading, temperatures, num_a_params);
        % pred = exp(isotherm_fun(params, loading, temperatures, num_a_params));
        residual = sum(weights .* (pressure - pred).^2);
        
    end
    function [fitted_params, RMSE, r2, exitflag] = run_ga(isotherm_opt, params)
        % Extracting data from isotherm options structure
        current_model = isotherm_opt.name;
        isotherm_fun = isotherm_opt.fun;
        params_lb = isotherm_opt.lb;
        params_ub = isotherm_opt.ub;
        num_a = isotherm_opt.num_a;
        num_params = isotherm_opt.num_params;
        
        % Fitness function handle
        fit_function = @(x)residual_func_ga(x, isotherm_fun, isotherm_data(:, 1), ...
                                        isotherm_data(:, 2), isotherm_data(:, 3), num_a, weight_vector);
        
        % options = optimoptions("ga", "Display", "off", "FunctionTolerance", 1e-14,...
        %                     "InitialPopulationMatrix", params, "MaxTime", 500);
        options = optimoptions("ga", "Display", "off", "FunctionTolerance", 1e-14,...
                                "MaxTime", 5000, "MaxGenerations", 1000000);
        try
            [fitted_params, residual, exitflag] = ga(fit_function, num_params, [], [], [], [], params_lb, params_ub,...
                 [], options);
           
            if exitflag == 0
                error("Maximum Generations for isotherm model %s reached! Terminating the solution...", ...
                        current_model);
            elseif exitflag < -1
                error("Solver failed to find the isotherm parameters for isotherm model %s",...
                        current_model);
            end

            % Root Mean Square Error
            residual_value = residual;  % Already summed errors
            
            RMSE = sqrt(residual_value./(N-isotherm_opt.num_params));

            % Coefficient of Determination
            ss_residual = residual;
            sq_var = sum((isotherm_data(:, 1)-mean(isotherm_data(:, 1), 1)).^2);
            r2 = 1 - ss_residual/sq_var;
        catch ME1
            disp(ME1.message);
            fitted_params = NaN(size(params_lb));
            RMSE = 1e05;
            r2 = 0;
            exitflag = -2;
        end
    end
%% Function for estimating covarience matrix and standard error
    function [cov_mat, std_error] = get_cov(fitted_parameters, isotherm_opt)
         % Extracting data from isotherm options structure
        isotherm_fun = isotherm_opt.fun;
        num_a = isotherm_opt.num_a;

        % Residual function handle
        fit_function = @(x)residual_func(x, isotherm_fun, isotherm_data(:, 1), ...
                                        isotherm_data(:, 2), isotherm_data(:, 3), num_a, weight_vector);
       
        options = optimset('Display', 'off', 'Diagnostics', 'off', 'MaxIter', 100000,...
                                'TolFun', 1e-12, 'TolX', 1e-12,...
                                'MaxFunEvals', 10000);
        try
            % Just a single run of lsqnonlin to get Jacobian
            [~,resnorm, ~,~,~,~,J] = lsqnonlin(fit_function, fitted_parameters, [], [], options);
        
            % Mean Squared Error
            MSE = resnorm/(N-isotherm_opt.num_params);
            
            % Covarience matrix
            cov_mat = MSE .* inv(J'*J);
            std_error = sqrt(diag(cov_mat));
        catch
            cov_mat = NaN * ones(size(isotherm_opt.num_params));
            std_error = NaN * ones(size(isotherm_opt.num_params), 1);
        end
    end
end

