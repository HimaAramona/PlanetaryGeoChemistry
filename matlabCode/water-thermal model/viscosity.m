function nu = viscosity(temp,f_water,g,rho_m,P)
%This function calculates the viscosity of the mantle following the
%parameterization of Sandu et al. (2011). The viscosity is dependent on the
%water abundance in the mantle, the temperature. The water abundance
%must be in mass fraction of the mantle.

%convert the mass fraction of the mantle into concentration in units of
%number of atoms of H per 10^6 Si atoms in olivine. Assume an olivine
%solid solution of (fayalite/forsterite) = 1/9
forsterite = 0.9;
fayalite = 0.1;
mfor = 140.0E-3;
mfay = 204.0E-3;
mH2O = 18.02e-3;
molv = forsterite * mfor + fayalite * mfay;
% rho_m = 3.3e3; %hold density at surface equal to uncompressed olivine

% P = rho_m * g * z; % pressure in Pa at 100 km (estimate for Db thickness)

%water abundance in OH molecules per 10^6 Si atoms
C_OH = f_water * molv * 1e6 * 2 / mH2O;

%constants relating C_OH to water fugacity (from Li et al. 2008)
c0 = -7.9859;
c1 = 4.3559;
c2 = -0.5742;
c3 = 0.0337;

%calculate ln f_H2O using the water abundance and the above parameters
% ln f_H2O = c0 + c1 lnCOH + c2 ln^2 COH + c3 ln^3 COH (Li et al. 2008)
logfH2O = c0 + c1 * log(C_OH) + c2 * (log(C_OH))^2 + c3 * (log(C_OH))^3;

%viscosity constants
    %calibration constant P-dependent
    eta_0 = 1.24e14;
%     eta_0 = 3e17; %P-independent
    %fugacity exponent
    r = 1.0;
    %activation energy for creep 
    Qa = 335e3; %J/mol
%     %activation volume
    V = 4e-6; %no pressure dependence of viscosity
%     V = 0.0;
    %ideal gas constant
    R = 8.31447; % J/mole/K
    %mantle viscosity

%effective dynamic viscosity (Pa s)
eta_eff = eta_0  * (exp(logfH2O))^-r * exp((Qa + P*V)/R/temp);
%kinematic viscosity (m^2 / s)
nu = eta_eff/rho_m; 
end