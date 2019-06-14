function [rsub,Dhyd] = regassing(qm,f_water,Rp,m0,Mm)
%This function calculates the regassing or subduction rate of volatiles
%into the mantle. It needs to solve first for the mass fraction of
%volatiles in the hydrated layer and the average hydrated layer thickness.

%the hydrated layer is serpentine, so the mass fraction of volatiles in the
%hydrated layer is the mass fraction of water in serpentine
f_h = 0.03;
chi_r = 0.03;
rho = 3.3e3;

%the hydrated layer is equal to the thickness of the boundary layer at
%which the temeprature reaches 700 C. The hydrated layer cannot be thicker
%than the boundary layer. (or can it?)
Ts = 300;
Tf = 700+273;
km = 4.2;

Dhyd = (Tf - Ts) * km / qm;

surface = (m0 - f_water * Mm);

if f_h*rho*4*pi()/3*(Rp^3 - (Rp-Dhyd)^3) > surface && surface > 0;
    Dhyd = Rp - ((Rp^3 - 3 * surface / f_h/rho/4/pi())^(1/3));
elseif surface <= 0
    Dhyd = 0;
end


rsub = f_h * Dhyd * rho * chi_r;

end




