classdef Isotherm
    properties
        isotherm_type
        Params
    end

    methods
        function obj = Isotherm(type_of_isotherm, parameters)
            % Object Constructor Method
            obj.isotherm_type = type_of_isotherm;
            obj.Params = parameters;
        end

        function loading = pure_loading(obj, P, T)
            % Returns the pure component loading values based on the Isotherm Model and Params.
            if strcmpi(obj.isotherm_type, 'DS-Langmuir')
                [~, loading] = obj.Langmuir_DS(P, T);

            elseif strcmpi(obj.isotherm_type, 'Langmuir-Freundlich')
                [~, loading] = obj.Langmuir_Freundlich(P, T);
            
            elseif strcmpi(obj.isotherm_type, 'BET')
                [~, loading] = obj.BET(P, T);

            elseif strcmpi(obj.isotherm_type, 'Quadratic')
                [~, loading] = obj.Quadratic(P, T);
            
            elseif strcmpi(obj.isotherm_type, 'Temkin')
                [~, loading] = obj.Temkin(P, T);
            
            else 
                error('Unidentified Isotherm Type: Check the Isotherm type in Isotherm');
            end
        end

        function reduced_grand_potential = reduced_grand_potential_calc(obj, P, T)
            % Returns the pure component reduced grand potential based on the Isotherm Model and Params.
            if strcmpi(obj.isotherm_type, 'DS-Langmuir')
                reduced_grand_potential = obj.Langmuir_DS(P, T);

            elseif strcmpi(obj.isotherm_type, 'Langmuir-Freundlich')
                reduced_grand_potential = obj.Langmuir_Freundlich(P, T);
            
            elseif strcmpi(obj.isotherm_type, 'BET')
                reduced_grand_potential = obj.BET(P, T);

            elseif strcmpi(obj.isotherm_type, 'Quadratic')
                reduced_grand_potential = obj.Quadratic(P, T);
            
            elseif strcmpi(obj.isotherm_type, 'Temkin')
                reduced_grand_potential = obj.Temkin(P, T);

            else 
                error('Unidentified Isotherm Type: Check the Isotherm type in Isotherm');
            end
        end
        
        %% FUNCTIONS FOR RESPECTIVE ISOTHERMS
        function [reduced_grand_potential, loading_Langmuir_DS] = Langmuir_DS(obj, P, T)
            % Site b
            m_b = obj.Params(1);
            b_0 = obj.Params(2);
            U_b = obj.Params(3);
            
            % Site d
            m_d = obj.Params(5);
            d_0 = obj.Params(6);
            U_d = obj.Params(7);
            R = 8.314          ;
            b = b_0 .* exp(-U_b/R./T);
            d = d_0 .* exp(-U_d/R./T); 

            reduced_grand_potential = m_b .* log(1 + b.*P) + m_d .* log(1 + d.*P);
            
            if nargout>1
                loading_Langmuir_DS = m_b .* b .*  P ./ (1 + b .* P) + m_d .* d .* P ./ (1 + d .* P);
            end
        end

        function [reduced_grand_potential, loading_Langmuir_Freundlich] = Langmuir_Freundlich(obj, P, T)
            % Site b
            m_b = obj.Params(1);
            b_0 = obj.Params(2);
            U_b = obj.Params(3);
            n_b = obj.Params(4);

            % Site d
            m_d = obj.Params(5);
            d_0 = obj.Params(6);
            U_d = obj.Params(7);
            n_d = obj.Params(8);

            R = 8.314          ;
            b = b_0 .* exp(-U_b/R./T);
            d = d_0 .* exp(-U_d/R./T);    
            
            reduced_grand_potential = m_b / n_b .* log(1 + b.*P.^n_b) + m_d / n_d .* log(1 + d.*P.^n_d);

            if nargout>1
                loading_Langmuir_Freundlich = m_b .* b .*  (P.^n_b) ./ (1 + b .* (P.^n_b)) + m_d .* d .* (P.^n_d) ./ (1 + d .* P.^n_d);
            end
        end

        function [reduced_grand_potential, loading_BET] = BET(obj, P, T)
            m = obj.Params(1);
            b_surface_0 = obj.Params(2);
            b_layers_0 = obj.Params(3);
            U_surface = obj.Params(4);
            U_layers = obj.Params(5);

            R = 8.314;
            b_surface = b_surface_0 .* exp(-U_surface/R./T);
            b_layers = b_layers_0 .* exp(-U_layers/R./T); 

            reduced_grand_potential = m .* log((1 + b_surface.*P - b_layers.*P) ./ (1 - b_layers.*P));

            if nargout>1
                loading_BET = m .* b_surface .* P ./ (1 - b_layers.*P) ./ (1 - b_layers.*P + b_surface.*P);
            end
        end

        function [reduced_grand_potential, loading_Temkin] = Temkin(obj, P, T)
            m = obj.Params(1);
            b = obj.Params(2);
            thetha = obj.Params(3);

            reduced_grand_potential = m .* (log(1 + b.*P) - 0.5*thetha .* (b.*P/(1 + b.*P)).^2);

            if nargout>1
                first_term = m .* (b.*P/(1 + b.*P));
                first_multiplier = m.*thetha .* (b.*P/(1 + b.*P)).^2;
                second_multiplier = (-1 + b.*P./(1 + b.*P));
                loading_Temkin = first_term + first_multiplier .* second_multiplier;
            end
        end

        function [reduced_grand_potential, loading_Quadratic] = Quadratic(obj, P, T)
            m = obj.Params(1);
            b = obj.Params(2);
            c = obj.Params(3);

            reduced_grand_potential = m .* log(1 + b.*P + c.*(P.^2));

            if nargout>1
                loading_Quadratic = m .* (b.*P + 2*c.*P.^2) ./ (1 + b.*P + c.*P.^2);
            end
        end
        
%         function reduced_grand_potential = reduced_grand_potential_calc(obj, pressure, temperature)
%             % Calculate the Reduced Grand Potential based on the Isotherm Model and Params
%             if strcmpi(obj.isotherm_type, 'DS-Langmuir')                
%                 % Site b
%                 m_b = obj.Params(1);
%                 b_0 = obj.Params(2);
%                 U_b = obj.Params(3);
%                 
%                 % Site d
%                 m_d = obj.Params(4);
%                 d_0 = obj.Params(5);
%                 U_d = obj.Params(6);
%                 R = 8.314          ;
%                 b = b_0 .* exp(-U_b/R./temperature);
%                 d = d_0 .* exp(-U_d/R./temperature); 
% 
%                 reduced_grand_potential = m_b .* log(1 + b.*pressure) + m_d .* log(1 + d.*pressure);  
%             
%             elseif strcmpi(obj.isotherm_type, 'Langmuir-Freundlich')
%                 % Site b
%                 m_b = obj.Params(1);
%                 b_0 = obj.Params(2);
%                 U_b = obj.Params(3);
%                 n_b = obj.Params(4);
% 
%                 % Site d
%                 m_d = obj.Params(5);
%                 d_0 = obj.Params(6);
%                 U_d = obj.Params(7);
%                 n_d = obj.Params(8);
% 
%                 R = 8.314          ;
%                 b = b_0 .* exp(-U_b/R./temperature);
%                 d = d_0 .* exp(-U_d/R./temperature);
%                 reduced_grand_potential = m_b / n_b .* log(1 + b .* pressure.^n_b) + m_d / n_d .* log(1 + d.*pressure.^n_d);
%             else
%                 error('Unidentified Isotherm Type: Check the Isotherm type in Isotherm');
%             end
%         end
    end
end