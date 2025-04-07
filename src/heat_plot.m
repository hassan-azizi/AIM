function [] = heat_plot(ax, pressure_data, loading_data, T_flag, T_array, isotherm_struc,...
                             xlogscale_flag, ylogscale_flag, unit_loading)    
    %% Function variables    
    line_width = 2.5;
    scatter_marker_size = 60;
    fontsize = 18;
    fontweight = 'bold';
    % fontname = 'AvantGrande';
    fontname = 'Calibri';    
    legend_array = {};      % Array to hold legend entries
    % color_array = ["blue", "green", "red"];       % Array to hold colors

    if ~isempty(isotherm_struc)
        %% Organizing data
        loading_data = rmmissing(reshape(loading_data, [], 1));
        % Removing 
        [idx_ret, ~] = find(loading_data~=0.0);
        loading_data = loading_data(idx_ret, 1);
        dH = isotherm_struc.dH ./ 1e03;  % Conversion from J/mol to kJ/mol
        % dH = -1.*dH;                     % For plotting purposes
        if isscalar(dH)
            dH = ones(size(loading_data)).*dH;
        elseif length(dH) ~= length(loading_data)
            error("Inconsistent size of loading data and enthalpy of adsorption.");
        end

        %% Plotting heat data
        cla (ax, "reset");
        
        % Plot positive dH with negative sign in ylabel 
        scatter(ax, loading_data, -1.0.*dH, scatter_marker_size, 'filled', 'o', 'MarkerEdgeColor','k');
           
        % Checking for log scale requirement
        if xlogscale_flag
            xscale(ax, "log");
        end
    
        if ylogscale_flag
            yscale(ax, "log");
        end
    
        ax.Box = "on";
        ax.Color = [1 1 1];
        
        % if ~isempty(unit_pressure)
        %     x_label = sprintf("Pressure (%s)", unit_pressure);
        % else
        %     x_label = sprintf("Pressure");
        % end
        % y_label = sprintf("Enthalpy of Adsorption dH (kJ/mol)");
        y_label = sprintf("Heat of Adsorption %sdH (kJ/mol)", char(8722));
    
        if ~isempty(unit_loading)
            x_label = sprintf("Gas Uptake (%s)", unit_loading);
        else
            x_label = sprintf("Gas Uptake");
        end
    
        xlabel(ax, x_label, FontSize = fontsize, FontName=fontname, FontWeight=fontweight);
        ylabel(ax, y_label, FontSize = fontsize, FontName=fontname, FontWeight=fontweight);
        
        % legend(ax, {'dH'}, Location="best");
        
        grid (ax, "on")
        ax.GridAlpha = 0.05;
        ytickformat(ax, "%.1f");

        ax.FontName = fontname;
        ax.FontWeight = fontweight;
        ax.FontSize = fontsize;
        % hold(ax, 'off');
    
        % % Ensuring consistent color scheme for plots
        % % Grab handles of scatter plots
        % scatter_handles = findobj(ax, 'Type', 'Scatter');
        % % if isempty(scatter_handles)
        % %     scatter_handles = 0;
        % % end
        % 
        % % Grab handles of line plots
        % line_handles = findobj(ax, 'Type', 'Line');
        % % if isempty(line_handles)
        % %     line_handles = 0;
        % % end
        % 
        % % Get color order from the axes
        % colors = ax.ColorOrder;
        % 
        % for i=1:length(scatter_handles)
        %     % Setting the corresponding scatter and line pair to have the same
        %     % color index
        %     if ~isempty(scatter_handles)
        %         set(scatter_handles(i), 'CData', colors(i, :));
        %     end
        % 
        %     if ~isempty(line_handles)
        %         set(line_handles(i), 'Color', colors(i, :));
        %     end
        % end
    end
end