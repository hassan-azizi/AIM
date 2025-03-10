function loading = EDSL_func(num_components, Pressure, gas_phase_mol_fraction, temperature, iso_params, T_flag)
    
    partial_pressures = Pressure .* gas_phase_mol_fraction;
    
    if T_flag
        T = temperature;
        R = 8.314;
        dH = iso_params(end-1, 1:num_components);
        T_ref = iso_params(end, 1:num_components);
        % Normalzing the partial pressure based on CC equation
        partial_pressures = partial_pressures .* exp(-1.0 * dH/R .* (1./T - 1./T_ref));             
    end
 
    b_sat_loading = iso_params(2, 1:num_components);
    b = iso_params(3, 1:num_components);
    d_sat_loading = iso_params(4, 1:num_components);
    d = iso_params(5, 1:num_components);

    site_b_denom = 1 + sum(b .* partial_pressures, 2);
    site_d_denom = 1 + sum(d .* partial_pressures, 2);
    loading_site_b = b_sat_loading .* b .* partial_pressures ./ site_b_denom;
    loading_site_d = d_sat_loading .* d .* partial_pressures ./ site_d_denom;
    
    loading_total = loading_site_b + loading_site_d;
    loading = loading_total;
end