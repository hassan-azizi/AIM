function loading = EDSL_func_mix_pred_app(num_components, iso_params, Pressure, gas_phase_mol_fraction, temperature)
    
    partial_pressures = Pressure .* gas_phase_mol_fraction;
    % T = temperature;
    % R = 8.314;
 
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