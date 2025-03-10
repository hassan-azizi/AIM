function isotherm_struc = Isotherm_fitting(isotherm_model, pressure_data, loading_data, guess_matrix, ...
                                           lb_matrix, ub_matrix, weight_assigment, plot_res,...
                                           algorithm_id, iso_models_available, P_saturation)
    %% Function Variables
    N = size(pressure_data, 1);     % number of data points
    
    % P_saturation is required for Klotz, Dubinin-Astakhov, and Do-Do model
    if nargin<11
        P_saturation = 1.0;
    end

    % Weights assignments correponding to the data points
    switch(weight_assigment)
        case 'Uniform'
            weight_vector = ones(size(pressure_data));
        case 'Biased'
            % Higher weight values for smaller loadings
            weight_vector = ones(size(pressure_data)) .* mean(loading_data) ./ loading_data;
        otherwise
            error("Please use the valid weight assignment method...")
    end
    
    % Checking if isotherm fitting is recommended for all isotherm
    % iso_models_available = {'SS-Langmuir', 'DS-Langmuir',...
    %                         'SS-Langmuir-Freundlich', 'DS-Langmuir-Freundlich',...
    %                          'Quadratic', 'Temkin', 'BET',...
    %                          'Sips', 'Toth','Dubinin-Astakhov', 'Klotz', 'Do-Do'};
    
    % In case of auto mode only the following general isotherm model will
    % be searched. S-Shaped isotherm can be chosen manually
    % iso_models_available = {'SS-Langmuir', 'DS-Langmuir',...
    %                             'SS-Langmuir-Freundlich', 'DS-Langmuir-Freundlich',...
    %                             'Quadratic', 'Temkin', 'BET',...
    %                             'Sips', 'Toth'};

    if strcmpi(isotherm_model, 'auto')
        isotherm_fitting_mode = 1;
        iso_model_array = iso_models_available;
    else
        isotherm_fitting_mode = 0;
        if isempty(isotherm_model)
            isotherm_model = 'SS-Langmuir';
        end
        iso_model_array = {isotherm_model};
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

    %% Isotherm Fitting
    if isotherm_fitting_mode
        RMSE_diff_iso = zeros(1, length(iso_model_array));
        r2_diff_iso = zeros(1, length(iso_model_array));
        exit_flag_diff_iso = ones(1, length(iso_model_array));
        fitted_param_struc_diff_iso = struct;

        for i = 1:length(iso_model_array) 
            current_isotherm_model = char(iso_model_array(i));
            
            % Calling Isotherm options
            isotherm_opt = isotherm_fit_opt(current_isotherm_model, loading_data);
            
            % Changing the default values of initial guesses/upper and
            % lower bounds to user specified values  
            idx = find(strcmpi(iso_models_available, current_isotherm_model));
            isotherm_opt.guess = guess_matrix(1:isotherm_opt.num_params, idx)';
            isotherm_opt.lb = lb_matrix(1:isotherm_opt.num_params, idx)';
            isotherm_opt.ub = ub_matrix(1:isotherm_opt.num_params, idx)';

            % Fitting using different algorithms
            [fitted_params_current_iso, min_RMSE_current_iso, r2_current_iso, exitflag_current_iso, ~] = fit_isotherm_diff_algo(isotherm_opt, fitting_algo_mode, fitting_algorithm_array);
            
            % Final best parameters and RMSE for current isotherm model
            RMSE_diff_iso(i) = min_RMSE_current_iso;
            r2_diff_iso(i) = r2_current_iso;
            exit_flag_diff_iso(i) = exitflag_current_iso;
            fitted_param_struc_diff_iso.(strrep(current_isotherm_model, '-', '_')) = fitted_params_current_iso;
        end
        
        % Choosing the best isotherm model
        [min_RMSE_iso, idxmin_iso] = min(RMSE_diff_iso);
        isotherm_min_RMSE = char(iso_model_array(idxmin_iso));
        fitted_params_best_iso = fitted_param_struc_diff_iso.(strrep(isotherm_min_RMSE, '-', '_'));
        min_RMSE_best_iso = min_RMSE_iso;
        r2_value_best_iso = r2_diff_iso(idxmin_iso);
        exitflag_best_iso = exit_flag_diff_iso(idxmin_iso);
        best_isotherm = isotherm_min_RMSE;
        
    else
        % Isotherm model already specified by user
        current_isotherm_model = char(iso_model_array);
        
        % Calling Isotherm options
        isotherm_opt = isotherm_fit_opt(current_isotherm_model, loading_data);
        
        % Changing the default values of initial guesses/upper and
        % lower bounds to user specified values  
        idx = find(strcmpi(iso_models_available, current_isotherm_model));
        isotherm_opt.guess = guess_matrix(1:isotherm_opt.num_params, idx)';
        isotherm_opt.lb = lb_matrix(1:isotherm_opt.num_params, idx)';
        isotherm_opt.ub = ub_matrix(1:isotherm_opt.num_params, idx)';
    
        % Fitting using different algorithms
        [fitted_params_current_iso, min_RMSE_current_iso, r2_current_iso, exitflag_current_iso, ~] = fit_isotherm_diff_algo(isotherm_opt, fitting_algo_mode, fitting_algorithm_array);
        
        fitted_params_best_iso = fitted_params_current_iso;
        min_RMSE_best_iso = min_RMSE_current_iso;
        r2_value_best_iso = r2_current_iso;
        exitflag_best_iso = exitflag_current_iso;
        best_isotherm = current_isotherm_model;      
    end
    
    %% Organizing output 
    if r2_value_best_iso <= 0 &&  exitflag_best_iso < 0
        error(['Failed to fit the isotherm. Try another isotherm model or ' ...
                'use better initial guesses...!']);
    end
    
    % Calling Isotherm options for the best isotherm model
    isotherm_opt = isotherm_fit_opt(best_isotherm, loading_data);
    
    % Get Covarience and standard error
    [~, std_error] = get_cov(fitted_params_best_iso, isotherm_opt);
    std_error = full(std_error);
    
    isotherm_struc.name = isotherm_opt.name;
    isotherm_struc.fun = isotherm_opt.fun;
    isotherm_struc.fitted_params = fitted_params_best_iso;
    isotherm_struc.RMSE = min_RMSE_best_iso;
    isotherm_struc.r2 = r2_value_best_iso;
    isotherm_struc.EF = exitflag_best_iso;
    isotherm_struc.T_flag = isotherm_opt.T_flag;
    isotherm_struc.p_sat = isotherm_opt.p_sat;
    isotherm_struc.std_error = std_error;

    if isotherm_fitting_mode
        isotherm_struc.RMSE_all = RMSE_diff_iso;
        isotherm_struc.r2_all = r2_diff_iso;
    end
    %
    %% Plotting
    if strcmpi(plot_res, "True")

        loading_pred = isotherm_struc.fun(isotherm_struc.fitted_params, pressure_data);
        
        scatter(pressure_data, loading_data);
        hold on
        
        plot(pressure_data, loading_pred)
        legend({'Data', best_isotherm})
        hold off
    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            %% MISCELLENEOUS FUNCTIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Residual Function
    function residual = residual_func(params, isotherm_fun, pressure, loading, weights)    
        pred = isotherm_fun(params, pressure);
        % residual = weights .* (loading - pred).^2;
         residual = (weights .* (loading - pred));
    end

    %% Isotherm Fitting Function
    function [fitted_params, resnorm, exitflag, RMSE, r2] = fit_isotherm(isotherm_opt, current_options, current_algo)
        
        % Extracting data from isotherm options structure
        current_model = isotherm_opt.name;
        isotherm_fun = isotherm_opt.fun;
        params_lb = isotherm_opt.lb;
        params_ub = isotherm_opt.ub;
        initial_guess = isotherm_opt.guess;
        p_sat_flag = isotherm_opt.p_sat;

        % Residual function handle
        if p_sat_flag
            % If p_sat_flag is required the pressure will be normalized
            % using P_saturation. Applicable to Klotz, Dubinin-Astakhov,
            % and Do-Do models.
            fit_function = @(x)residual_func(x, isotherm_fun, pressure_data./P_saturation, loading_data, weight_vector);
        else
            fit_function = @(x)residual_func(x, isotherm_fun, pressure_data, loading_data, weight_vector);
        end

        try
            [fitted_params, resnorm, ~, exitflag] = lsqnonlin(fit_function, initial_guess, params_lb, params_ub, current_options);
        
            if exitflag == 0
                error("Maximum Iterations for isotherm model %s and algorithm %s reached! Terminating the solution...", ...
                        current_model, current_algo);
            elseif exitflag < 0
                error("Solver failed to find the isotherm parameters for isotherm model %s and algorithm %s",...
                        current_model, current_algo);
            end

            % Root Mean Square Error
            % residual_value = sum(fit_function(fitted_params));
            residual_value = resnorm; 
            RMSE = sqrt(residual_value./(N-isotherm_opt.num_params));

            % Coefficient of Determination
            ss_residual = resnorm;
            % ss_residual = sum(fit_function(fitted_params).^2);
            sq_var = sum((loading_data-mean(loading_data, 1)).^2);
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
                options = optimset('Display', 'off', 'Diagnostics', 'off', 'MaxIter', 100000,...
                                'Algorithm', current_algo, 'TolFun', 1e-12,...
                                'TolX', 1e-12, 'MaxFunEvals', 10000);
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

        % Fitting using user-specified algorithm
        else
            current_algo = char(fitting_algorithm_array);
            % Fitting options
            options = optimset('Display', 'off', 'Diagnostics', 'off', 'MaxIter', 100000,...
                                'Algorithm', current_algo, 'TolFun', 1e-12, 'TolX', 1e-12,...
                                'MaxFunEvals', 10000);
            [fitted_params, ~, exitflag, min_RMSE, r2_value] = fit_isotherm(isotherm_opt, options, current_algo);
            fitted_params_current_iso = fitted_params;
            min_RMSE_current_iso = min_RMSE;
            exitflag_current_iso = exitflag;
            r2_value_current_iso = r2_value;
            best_algorithm = current_algo;
        end
    end
    %% Function for estimating covarience matrix and standard error
    function [cov_mat, std_error] = get_cov(fitted_parameters, isotherm_opt)
        % Extracting data from isotherm options structure
        isotherm_fun = isotherm_opt.fun;
        % params_lb = isotherm_opt.lb;
        % params_ub = isotherm_opt.ub;
        p_sat_flag = isotherm_opt.p_sat;

        % Residual function handle
        if p_sat_flag
            % If p_sat_flag is required the pressure will be normalized
            % using P_saturation. Applicable to Klotz, Dubinin-Astakhov,
            % and Do-Do models.
            fit_function = @(x)residual_func(x, isotherm_fun, pressure_data./P_saturation, loading_data, weight_vector);
        else
            fit_function = @(x)residual_func(x, isotherm_fun, pressure_data, loading_data, weight_vector);
        end
        
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
        catch ME1
            cov_mat = NaN * ones(isotherm_opt.num_params, isotherm_opt.num_params);
            std_error = NaN * ones(isotherm_opt.num_params, 1);
        end
    end
end