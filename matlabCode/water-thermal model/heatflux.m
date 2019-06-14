function [qm,Db,uc,Ra,nu] = heatflux(temp,f_water,Rp,Rc,g,avgfact,rho_m,avgP)
%This function calculates the mantle heat flux, which depends on the
%temperature and viscosity of the mantle. The viscosity depends on
%temperature, pressure, and water abundance and is calculated with the
%function viscosity.m

%mantle thermal conductivity (W/m/K)
km = 4.2;
%depth of the convective zone
Z = Rp - Rc;
%critical Rayleigh number
Ra_cr = 1.1e3;
%exponent
beta = 0.33;
% coefficient of thermal expansion (1/K)
alpha = 2e-5;
%mantle heat capacity (J/kg/K)
cp = 1.2e3; 
%mantle thermal diffusivity (m^2 / s)
kappa = km / rho_m / cp;
%surface temperature (K)
Ts = 300;
Tp = temp/avgfact;

%find the viscosity (m^2/s)
nu = viscosity(temp,f_water,g,rho_m,avgP);

%Calculate Rayleigh number
Ra = (g * alpha * (temp - Ts) * Z^3) / nu / kappa;

%boundary layer depth (m)
Db = Z * (Ra_cr / Ra)^beta;

%find temperature at the base of the boundary layer
Tb = Tp + Tp * (alpha * g * Db / cp);

%calculate the heat flux (W/m^2)
qm = km * (Tb - Ts) / Db;

%spreading time (s)
ts = Db^2 / (5.38 * kappa);%s
%spreading velocity (m/s)
uc = Z/ts; %m/s

end

