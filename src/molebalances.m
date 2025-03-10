function molar_accumulation = molebalances(P, T, y, Ph, Th, yh, vh,...
                                dydz, d2ydz2, dpdz, dTdz, dPdt, dTdt, dxdt,...
                                 epsilon, epsilon_3, Pe, phi, dz, N)
%% This function returns the molar accumulation of chemical specie
%   Input:
%   P = Dimensionless Pressure
%   T = Dimensionless Temperature
%   y = Mole fraction vector 
%   Ph = Dimensionless Pressure at walls of node
%   Th = Dimensionless Temperature at walls of node
%   yh = Molefraction at walls of node
%   vh = Velocity at walls of node
%   r_sp = Net Rate of reaction of specie
%   dydz, d2ydz2, dpdz, dTdz = Spatial gradients of mole fraction, pressure
%                               and temperature respectively.
%   dPdt, dTdt, dx1dt = Temporal gradients of pressure, temperature and 
%                           adsorbion.
%   epsilon = Bed Porosity
%   epsilon_3 = Total porosity
%   Pe = Peclet number
%   phi = Adsorbption paramter
%   phi_r = Reaction parameters
%   dz = differential length of the column
%   N = Nodes Index

%   Output:
%   rate of change of molar fraction in node
                                
    %%  
%   *4) Component Mass Balance (Based on Mole Fraction)*
%   
    dydt1 = zeros(size(y))   ;
    dydt2 = zeros(size(y))   ;
    dydt3 = zeros(size(y))   ;
    dydt4 = zeros(size(y))   ;
    dydt5 = zeros(size(y))   ;
%%  
%   $$ \frac{\partial y}{\partial \tau} = \frac{1}{Pe} \big(\frac{{
%      \partial}^2 y}{\partial {Z}^2}+\frac{1}{\bar{P}}\frac{\partial
%      \bar{P}}{\partial Z}\frac{\partial y}{\partial Z}-\frac{1}{\bar{T}}
%      \frac{\partial \bar{T}}{\partial Z}\frac{\partial y}{\partial Z}
%      \big)-\bar{v}\frac{\partial y}{\partial Z}+\frac{\Psi \bar{T}}{
%      \bar{P}} \big((y-1)\frac{\partial \bar{x_{1}}}{\partial \tau}+y
%      \frac{\partial\bar{x_{2}}}{\partial \tau}\big) $$
%   
%%  
%   4.1) Calculate the change in mole fraction due to diffusion
%   
%%  
%   $$ \frac{1}{Pe} \big(\frac{{\partial}^2 y}{\partial {Z}^2}+\frac{1}{
%      \bar{P}}\frac{\partial \bar{P}}{\partial Z}\frac{\partial y}{
%      \partial Z}-\frac{1}{\bar{T}}\frac{\partial \bar{T}}{\partial Z}
%      \frac{\partial y}{\partial Z} \big) $$
%   
%%  
    dydt1(2:N+1) = (epsilon / epsilon_3) * (1/Pe) .* (d2ydz2(2:N+1) +...
                    (dydz(2:N+1) .* dpdz(2:N+1)./P(2:N+1)) - (dydz(2:N+1)... 
                    .* dTdz(2:N+1)./ T(2:N+1)))                        ;
%   
%%  
%   4.2) Calculate the change in mole fraction due to advection
%   
%%  
%   $$ -\bar{v}\frac{\partial y}{\partial Z} $$
%   
%%  
    ypvt         = yh(1:N+1).*Ph(1:N+1).*vh(1:N+1)./Th(1:N+1)        ;
    PvT          = Ph(1:N+1).*vh(1:N+1)./Th(1:N+1)        ;       
    dydt2(2:N+1) = -(epsilon / epsilon_3) .* (T(2:N+1)./P(2:N+1))...
                        .*((ypvt(2:N+1)-ypvt(1:N)))./dz    ;
%     dydt2(2:N+1) = -(T(2:N+1)./P(2:N+1)).*((ypvt(2:N+1)-ypvt(1:N))... 
%                   -y(2:N+1).*(PvT(2:N+1)-PvT(1:N)))./dz             
    
%%  
%   4.3) Calculate the change in mole fraction due to adsorption/desorption  
    dydt3(2:N+1) = -(phi*T(2:N+1)./P(2:N+1)) .* dxdt(2:N+1);
%     dydt3(2:N+1) = (phi*T(2:N+1)./P(2:N+1)) .* ((y(2:N+1) - 1) .* dxdt(2:N+1)...
%                       + y(2:N+1) .* res_dxdt(2:N+1));                                     
%%
%   4.5) Mole fraction change due to Pressure and Temperature change
    dydt5(2:N+1) = -(y(2:N+1) ./ P(2:N+1)) .* dPdt(2:N+1) + (y(2:N+1) ./ T(2:N+1))...
                .* dTdt(2:N+1) ;
%% 
%   4.6) Total sum of all mole fraction changes
    dydt(2:N+1) = dydt1(2:N+1) + dydt2(2:N+1) + dydt3(2:N+1) + dydt4(2:N+1)...
                    + dydt5(2:N+1);


    molar_accumulation = dydt(2:N+1) ;

end