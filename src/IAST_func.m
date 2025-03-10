function partial_loadings = IAST_func(num_components, isotherm_params_array, Pressure, gas_phase_mol_fraction, T_array, T_flag, initial_guess)
    %% Parameter Unpacking and Decalarations
    size_of_pressure_vector = length(Pressure);     % Length of mole fraction vector
    partial_pressures = Pressure .* gas_phase_mol_fraction;
    adsorbed_mole_fractions = zeros(size_of_pressure_vector, num_components);
    %   
    %% Initialization of Guess for IAST Equations Solution
    % Solve the non-linear equations returened by the 'func' functions using fsolve method.
    % If Initial Guess has not been provided in input arguments, 
    % it will be guessed proportional to pure compoenent loading.
    
    if nargin<7 || ~(any(sum(initial_guess, 2)))
        loading_array = Isotherm_functions(num_components, isotherm_params_array, partial_pressures, T_array, T_flag, 1);
        loading_array(loading_array == 0) = 1e-07;
        initial_guess = loading_array ./ (sum(loading_array, 2));
    end
    
    % Use to avoid undefined values when the column is filled 
    % with only non-adsorbing gas or when the gas loading is very low
    tolerance = 1e-8;
    
    guess = initial_guess(:, 1:num_components-1);

%     options = optimset('Display','off', 'TolFun', 1e-10, 'MaxIter', 400, 'JacobPattern', jp);   % Options struc for fsolve
%     options = optimset('Display','off', 'TolFun', 1e-5, 'MaxIter', 400);   % Options struc for fsolve
    options = optimset('TolFun', 1e-06, 'Display', 'off', 'Diagnostics', 'off', 'TolX', 1e-09, 'FinDiffRelStep', 1e-9, 'Algorithm','levenberg-marquardt',...
        'MaxIter', 10e5);
    
    for k = 1:size_of_pressure_vector
        try
            % In case of noadsorbing component
            if sum(gas_phase_mol_fraction(k, :), 2) <= tolerance
                adsorbed_mole_fractions(k, :) = zeros(1, num_components)+tolerance;
            
            % In case of single adsorbing component, the initial guess will give the correct loading
            elseif num_components == 1
                adsorbed_mole_fractions(k, :) = initial_guess(k, :);
            
            % Calculate the IAST loadings 
            else               
                function_handle = @(x)residual_func(x, isotherm_params_array, partial_pressures(k, :), T_array(k, 1));
                [solution, ~, exitflag] = fsolve(function_handle, guess(k, :), options);
                % [solution, ~, ~, exitflag] = lsqnonlin(function_handle, guess(k, :), ...
                %                                       zeros(1, 1:num_components-1), ones(1, 1:num_components-1), options);
%                 [solution, ~, exitflag] = fsolve(function_handle, guess, options);

                if exitflag == 0
                    error("Maximum Iterations for fsolve reached! Terminating the solution.");
                elseif exitflag < 0
                    error("fsolve failed to solve the equations!")
                end
                adsorbed_mole_fractions(k, :) = [solution, 1 - sum(solution, 2)];

                ads_mol_frac_unity_tol = 1e-10;    % Tolerance value to check that sum of adsorbed mole fraction is unity
                if ~(ismembertol(sum(adsorbed_mole_fractions(k, :)), 1, ads_mol_frac_unity_tol))
                    error("The sum of adsorbed mole fractions is not unity, fsolve failed to solve the equations!")
                end
            end
        catch ME1
                error(['Failed to solve the IAST equations...!', ...
                       'Error:%s'], ME1.message);    
        end
    end
    
    adsorbed_mole_fractions(adsorbed_mole_fractions == 0) = tolerance;
    pressure_0 = partial_pressures ./ adsorbed_mole_fractions;
    loading_array = Isotherm_functions(num_components, isotherm_params_array, pressure_0, T_array, T_flag, 1);
    inverse_loading = sum(adsorbed_mole_fractions./loading_array, 2);
    loading_total = 1 ./ inverse_loading;
    partial_loadings = (adsorbed_mole_fractions .* loading_total);

    function residual = residual_func(adsorbed_MF, isotherm_params_array, partial_pressure, T)
        % Calculate the difference in spreading pressure for (N-1) components.
        % The spreading pressures are all calculated at the fictious pressure calculated based on the
        % adsorbed mole fraction and partial pressure of the components.
        % Fictitous Pressure = p_* = p_i / x_i
        
        s = 1e-08;       
        x = [(adsorbed_MF), (1-sum(adsorbed_MF))];
        x(x==0) = s;
        
        spreading_pressures = Isotherm_functions(num_components, isotherm_params_array, partial_pressure./x, T, T_flag, 0);
        
        spreading_pressures_diff = spreading_pressures(:, 1:num_components-1) - spreading_pressures(:, 2:num_components);
        if ~isreal(spreading_pressures_diff)
            residual = 1e05;
        else
            residual = spreading_pressures_diff;
        end
        if any(isnan(residual))
            error("NaN values in IAST: residual function...")
        end
    end
end