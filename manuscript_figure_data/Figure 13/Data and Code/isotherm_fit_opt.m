function opt = isotherm_fit_opt(isotherm_model, loading_data, Pressure)
    %% Function variables
    num_param_vector = [2, 4, 3, 6, 3, 3, 3, 3, 3, 6, 3, 4, 6];
    UL_general = 1e06;   % it was 1e03 initially, changed after ding STA model
    UL_langmuir_constant = 1;

    LL_theta = -100;
    UL_theta = 100;
    SCALE = 1.0;

    sat_loading_guess = SCALE * max(loading_data, [], "all");
    langmuir_k_guess = 1e-10;
    exponent_guess = 1.0;
    theta_guess = 0.0;
    BET_tol = 1e-8;       % Need tolerance to avoid getting undefined values in BET  
    
    if nargin<3
        P_tr_guess = 1.0;
    else
        try
            % Add code to find inflection point pressure
            P_tr_guess = get_P_inflection(Pressure, loading_data);
        catch
            P_tr_guess = mean(Pressure);
        end

        if isnan(P_tr_guess)
            P_tr_guess = mean(Pressure);
        end
    end
    %% Isotherm options structure based on the chosen isotherm model
    switch(isotherm_model)
        case 'SS-Langmuir'
            % opt.fun = @(params, P) (params(1)*params(2).*P) ./ (1+params(2).*P);
            opt.fun = @SS_Langmuir;
            opt.num_params = num_param_vector(1); 
            % Lower Bounds
            opt.lb = zeros(1, opt.num_params);
            % Upper Bounds
            opt.ub = UL_general .* ones(1, opt.num_params);
            opt.ub(2) = UL_langmuir_constant;       % Upper bound for langmuir constant b 
            opt.guess = [sat_loading_guess, langmuir_k_guess];
            opt.T_flag = 1;
            opt.p_sat = 0;

        case 'DS-Langmuir'
            % opt.fun = @(params, P) (params(1)*params(2).*P) ./ (1+params(2).*P) + (params(3)*params(4).*P) ./ (1+params(4).*P);
            opt.fun = @DS_Langmuir;
            opt.num_params = num_param_vector(2);
            % Lower Bounds
            opt.lb = zeros(1, opt.num_params); 
            % Upper Bounds
            opt.ub = UL_general .* ones(1, opt.num_params);
            opt.ub([2, 4]) = UL_langmuir_constant;  % Upper bound for langmuir constant b and d
            opt.guess = [0.90*sat_loading_guess, langmuir_k_guess,... 
                         0.10*sat_loading_guess, 0.5*langmuir_k_guess];
            opt.T_flag = 1;
            opt.p_sat = 0;

        case 'SS-Langmuir-Freundlich'
            % opt.fun = @(params, P) (params(1)*params(2).*P.^params(3)) ./ (1+params(2).*P.^params(3));
            opt.fun = @SS_Langmuir_Freundlich;
            opt.num_params = num_param_vector(3);
            % Lower Bounds
            opt.lb = zeros(1, opt.num_params);
            % Upper Bounds
            opt.ub = UL_general .* ones(1, opt.num_params);
            opt.ub(2) = UL_langmuir_constant;       % Upper bound for langmuir constant b
            opt.guess = [sat_loading_guess, langmuir_k_guess, exponent_guess];
            opt.T_flag = 1;
            opt.p_sat = 0;

        case 'DS-Langmuir-Freundlich'
            % opt.fun = @(params, P) (params(1)*params(2).*P.^params(3)) ./ (1+params(2).*P.^params(3)) ...
            %                     + (params(4)*params(5).*P.^params(6)) ./ (1+params(5).*P.^params(6));
            opt.fun = @DS_Langmuir_Freundlich;
            opt.num_params = num_param_vector(4);
            % Lower Bounds
            opt.lb = zeros(1, opt.num_params);
            % Upper Bounds
            opt.ub = UL_general .* ones(1, opt.num_params);
            opt.ub([2, 5]) = UL_langmuir_constant;  % Upper bound for langmuir constant b
            opt.guess = [0.75*sat_loading_guess, langmuir_k_guess, exponent_guess ... 
                         0.25*sat_loading_guess, 0.5*langmuir_k_guess, exponent_guess];
            opt.T_flag = 1;
            opt.p_sat = 0;

        case 'Quadratic'
            % opt.fun = @(params, P) params(1) .* (params(2).*P + 2*params(3).*P.^2) ./ (1 + params(2).*P + params(3).*P.^2);
            opt.fun = @Quadratic;
            opt.num_params = num_param_vector(5);
            % Lower Bounds
            opt.lb = zeros(1, opt.num_params);
            % Upper Bounds
            opt.ub = UL_general .* ones(1, opt.num_params);
            opt.ub([2, 3]) = UL_langmuir_constant;
            opt.guess = [0.5*sat_loading_guess, langmuir_k_guess, langmuir_k_guess^2];
            opt.T_flag = 0;
            opt.p_sat = 0;

        case 'Temkin'
            % opt.fun = @(params, P) (params(1)*params(2).*P) ./ (1+params(2).*P) ...
            %                + (params(1)*params(3) .* ((params(2).*P) ./ (1+params(2).*P)).^2) ...
            %                .* ((params(2).*P) ./ (1+params(2).*P) - 1);
            opt.fun = @Temkin;
            opt.num_params = num_param_vector(6);
            % Lower Bounds
            opt.lb = zeros(1, opt.num_params); 
            opt.lb(3) = LL_theta; % Because theta parameter can be negative
            % Upper Bounds
            opt.ub = UL_general .* ones(1, opt.num_params);
            opt.ub(2) = UL_langmuir_constant;
            opt.ub(3) = UL_theta;
            opt.guess = [sat_loading_guess, langmuir_k_guess, theta_guess];
            opt.T_flag = 0;
            opt.p_sat = 0;

        case 'BET'
            % opt.fun = @(params, P) (params(1)*params(2).*P) ./ (1 - params(3).*P) ./ (1 - params(3).*P + params(2).*P);
            opt.fun = @BET;
            opt.num_params = num_param_vector(7);
            opt.lb = zeros(1, opt.num_params);
            opt.ub = UL_general .* ones(1, opt.num_params);
            opt.ub([2, 3]) = UL_langmuir_constant .* BET_tol;
            opt.guess = [sat_loading_guess, langmuir_k_guess, langmuir_k_guess^2];
            % opt.guess = [sat_loading_guess, sat_loading_guess, langmuir_k_guess^2];
            opt.T_flag = 0;
            opt.p_sat = 0;

        case 'Sips'
            opt.fun = @Sips;
            opt.num_params = num_param_vector(8);
            opt.lb = zeros(1, opt.num_params);
            opt.ub = UL_general .* ones(1, opt.num_params);
            opt.ub(2) = UL_langmuir_constant;
            opt.guess = [sat_loading_guess, langmuir_k_guess, exponent_guess];
            % opt.guess = [sat_loading_guess, sat_loading_guess, langmuir_k_guess^2];
            opt.T_flag = 0;
            opt.p_sat = 0;

         case 'Toth'
            opt.fun = @Toth;
            opt.num_params = num_param_vector(9);
            opt.lb = zeros(1, opt.num_params);
            opt.ub = UL_general .* ones(1, opt.num_params);
            opt.ub(2) = UL_langmuir_constant;
            opt.guess = [sat_loading_guess, langmuir_k_guess, exponent_guess];
            % opt.guess = [sat_loading_guess, sat_loading_guess, langmuir_k_guess^2];
            opt.T_flag = 0;
            opt.p_sat = 0;

         case 'Structural-Transition-Adsorption'
            opt.fun = @STA;
            opt.num_params = num_param_vector(10);
            opt.lb = zeros(1, opt.num_params);
            opt.ub = UL_general .* ones(1, opt.num_params);
            % opt.ub([2, 3]) = UL_langmuir_constant .* BET_tol;
            % opt.guess = [sat_loading_guess, K_guess, C_guess, n_guess];
            opt.guess = [sat_loading_guess, langmuir_k_guess,...
                         sat_loading_guess, langmuir_k_guess,...
                         exponent_guess, P_tr_guess];
            opt.T_flag = 0;
            opt.p_sat = 0;
        
        case 'Dubinin-Astakhov'
            opt.fun = @DA_Isotherm;
            opt.num_params = num_param_vector(11);
            opt.lb = zeros(1, opt.num_params);
            opt.lb(2) = 1e-05;          % TO avoid inifnitny 
            opt.ub = UL_general .* ones(1, opt.num_params);
            % opt.ub([2, 3]) = UL_langmuir_constant .* BET_tol;
            opt.guess = [sat_loading_guess, 1, exponent_guess];
            opt.T_flag = 0;
            opt.p_sat = 1;

        case 'Klotz'
            opt.fun = @Klotz;
            opt.num_params = num_param_vector(12);
            opt.lb = zeros(1, opt.num_params);
            opt.ub = UL_general .* ones(1, opt.num_params);
            % opt.ub([2, 3]) = UL_langmuir_constant .* BET_tol;
            opt.guess = [sat_loading_guess, 1, 1, exponent_guess];
            opt.T_flag = 0;
            opt.p_sat = 1;

        case 'Do-Do'
            opt.fun = @DD_Isotherm;
            opt.num_params = num_param_vector(13);
            opt.lb = zeros(1, opt.num_params);
            opt.ub = UL_general .* ones(1, opt.num_params);
            % opt.ub([2, 3]) = UL_langmuir_constant .* BET_tol;
            % opt.guess = [sat_loading_guess, K_guess, C_guess, n_guess];
            opt.guess = [sat_loading_guess, 0.0, 0.0, 0.0,...
                         exponent_guess, 2*exponent_guess];
            opt.T_flag = 0;
            opt.p_sat = 1;
            
        case 'Virial'
            opt.fun = @Virial_func;
            opt.num_a = 7;      % Refers to number of "a" parameters in Virial equation
            opt.num_b = 3;      % Refers to number of "b" parameters in Virial equation
            opt.num_params = opt.num_a + opt.num_b;
            % opt.lb = zeros(1, opt.num_params);
            opt.lb = -Inf.* ones(1, opt.num_params);
            opt.ub = Inf .* ones(1, opt.num_params);
            % opt.ub([2, 3]) = UL_langmuir_constant .* BET_tol;
            opt.guess = zeros(1, opt.num_params);
            opt.T_flag = 0;
            opt.p_sat = 0;

        otherwise
            error('Invalid isotherm model: %s', isotherm_model);
    end

    opt.name = isotherm_model;

    %% Checking the initial guess values
    % if nargin > 2  % Initial guess specified by user
    %     %% Validation of initial guess
    %     % Checking if number of initial guesses/lower bounds/upper bounds is equal to number of isotherm parameter
    %     if size(initial_guess, 2) ~= opt.num_params
    %         warning("Number of initial guesses/lower bounds/upper bounds/ does not match the number of parameters for the chosen isotherm model. " + ...
    %                 "Using default initialization...")
    %         % opt.guess = default_initial_guess .* ones(1, opt.num_params);
    % 
    %     % Specifying the new bounds and checking if any initial guess value is not within bounds
    %     else
    %         % opt.guess = initial_guess;
    %         is_in_range = initial_guess >= opt.lb & initial_guess <= opt.ub;
    %         if ~all(is_in_range)
    %             warning("Some of the initial guess values provided are " + ...
    %                 "not within the bounds. Those values will be initialized using default initial guesses.");
    % 
    %             % index_out_range = find(~is_in_range);
    %             % opt.guess(index_out_range) = default_initial_guess .* ones(1, length(index_out_range));
    %         end
    % 
    %         index_in_range = find(is_in_range);
    %         opt.guess(index_in_range) = initial_guess(is_in_range);
    %     end      
    % end
    %
    %% Miscelleneous Functions
    % SS-Langmuir
    function loading = SS_Langmuir(params, P)
        loading = (params(1)*params(2).*P) ./ (1+params(2).*P);
    end
    
    % DS-Langmuir
    function loading = DS_Langmuir(params, P)
        loading = (params(1)*params(2).*P) ./ (1+params(2).*P)... 
                + (params(3)*params(4).*P) ./ (1+params(4).*P);
    end

    % SS-Langmuir-Freundlich
    function loading = SS_Langmuir_Freundlich(params, P)
        loading = (params(1)*params(2).*P.^params(3)) ./ (1+params(2).*P.^params(3));
    end

    % DS-Langmuir-Freundlich
    function loading = DS_Langmuir_Freundlich(params, P)
        loading = (params(1)*params(2).*P.^params(3)) ./ (1+params(2).*P.^params(3))...
                + (params(4)*params(5).*P.^params(6)) ./ (1+params(5).*P.^params(6));
    end

    % Quadratic Isotherm
    function loading = Quadratic(params, P)
        loading = params(1) .* (params(2).*P + 2*params(3).*P.^2) ...
                ./ (1 + params(2).*P + params(3).*P.^2);
    end

    % Temkin Isotherm
    function loading = Temkin(params, P)
        m = params(1);
        b = params(2);
        theta = params(3);

        lang_term = b.*P ./ (1+b.*P); 
        first_term = m .* lang_term;
        first_multiplier = m.*theta .* (lang_term).^2;
        second_multiplier = lang_term - 1;
        loading = first_term + first_multiplier .* second_multiplier;
    end
    
    % BET Isotherm
    function loading = BET(params, P)
        m = params(1);
        b_surface = params(2);
        b_layers = params(3);

        loading = m.*b_surface.*P ./ (1 - b_layers.*P) ./ (1 - b_layers.*P + b_surface.*P);
    end

    % Sips Isotherm
    function loading = Sips(params, P)
        m_sat = params(1);
        b = params(2);
        n = params(3);
        num = m_sat .* (b.*P).^(1/n);
        den = 1 + (b.*P).^(1/n);
        loading = num ./ den;
    end
    
    % Toth Isotherm
    function loading = Toth(params, P)
        m_sat = params(1);
        b = params(2);
        n = params(3);
        num = m_sat .* b.*P;
        den = (1 + (b.*P).^(n)).^(1/n);
        loading = num ./ den;
    end

    % Klotz Isotherm
    function loading = Klotz(params, P)
        % Here P = relative pressure
        m_sat = params(1);
        K = params(2);
        C = params(3);
        n = params(4);

        q = K.*P;
    
        num = m_sat.*C.*q .* (1 - (1+n).*q.^n + n.*q.^(n+1));
        den = (1-q) .* (1 + (C-1).*q - C.*q.^(n+1));        
        loading = num ./ den;
    end
    
    % DD-Isotherm
    function loading = DD_Isotherm(params, P)
        % Here P = relative pressure
        m_sat = params(1);
        f = params(2);
        K_f = params(3);
        K_mu = params(4);
        alpha = params(5);
        beta = params(6);
        
        num = f .* (K_f.*P) .* (1 - (1+beta).*P.^beta + beta.*P.^(beta+1));
        den = (1-P) .* (1 + (K_f-1).*P - K_f.*P.^(beta+1));

        first_term = num ./ den;
        second_term = (1-f) .* (K_mu.*P.^alpha) ./ (1 + K_mu.*P.^alpha);

        loading = m_sat.*(first_term + second_term);
    end

    % DA isotherm
    function loading = DA_Isotherm(params, P)
        % Here P = relative pressure
        m = params(1);
        epsilon = params(2);
        n = params(3);

        loading = m .* exp(-(1/epsilon .* log(1./P)).^n);
    end

    % STA Isotherm
    function loading = STA(params, P)
        q_NP = params(1);
        b_NP = params(2);
        q_LP = params(3);
        b_LP = params(4);
        s = params(5);
        P_tr = params(6);
        
        den_NP = (1+b_NP.*P);
        den_LP = (1+b_LP.*P);

        y = ((1+b_NP*P_tr)./(den_NP)).^q_NP .* ((1+b_LP*P_tr)./(den_LP)).^q_LP;

        sigma = y.^s ./ (1+y.^s);
        
        lang_NP = q_NP.*b_NP.*P./den_NP;
        lang_LP = q_LP.*b_LP.*P./den_LP;

        loading = (1-sigma).*lang_NP + sigma.*lang_LP;
    end
    
    % Virial polynomial function
    function pressure = Virial_func(params, loading, temperature, num_a_params)
        % returns log of pressure values
        % Flip function is required because the polyval accepts the
        % coefficients in the order of decreasing degree.
        a_terms = 1./temperature .* polyval(flip(params(1:num_a_params)), loading);
        b_terms = polyval(flip(params(num_a_params+1:end)), loading);
        pressure = a_terms + b_terms + log(loading);
    end

     
    %% Function to find the inflection point for STA isotherm
    function P_inflection = get_P_inflection(Pressure, loading)
        % smoothing the Pressure values and using spline interpolation
        size_p = length(Pressure);
        resolve_factor = 2;
        pressure_smoothed = linspace(min(Pressure), max(Pressure), resolve_factor*size_p);
        
        % Smoothing the loading using 'loess' model
        loading = smooth(loading, 'loess');
        loading_smoothed = spline(Pressure, loading, pressure_smoothed);
        
        % Derivative Calculation using Finite Differences
        dqdP = diff(loading_smoothed) ./ diff(pressure_smoothed);
        d2qdp2 = diff(dqdP) ./ diff(pressure_smoothed(2:end));
        
        % Idx where decond derivative changes sign, in case of multiple
        % inflection points, only the first one is selected
        inflection_idx = find(diff(sign(d2qdp2)), 1, "first");
        
        P_inflection = pressure_smoothed(inflection_idx+2); % +2 because of finite differences
    end
end