function deriv = derivatives(t,T_water,Mm,Rp,Rc,g,avgfact,m0,avgP)
%This function calculates the temperature derivative for the thermal
%evolution model deriv(1) and the rate of water exchange between the
%surface and the mantle deriv(2). This will be used with an ODE solver. The
%variables to solve for are T_water(1) = mantle temperature and T_water(2)
%= water abundance.

deriv = zeros(2,1);
t

%mantle density (kg/m3)
rho_m = 3 * Mm / (4 * pi() * (Rp^3 - Rc^3));
%heat capacity(J/kg/K)
cp = 1.2e3; 
% mid-ocean ridge length = 1.5x planet circumference (m)
Lridge = 2 * pi() * Rp * 1.5;
%Mass of the mantle (kg)
% Mmantle = 4.06e24; 
delta = m0 - T_water(2);
deriv2est = delta / 1e7; %estimated maximum water derivative
if T_water(2) < m0
    f_water = T_water(2) / Mm;
elseif T_water(2) >= m0
    T_water(2) = m0;
    f_water = m0/Mm;
end

% Radioactive heat constants (W/kg), from Schubert et al. (2001) Ch. 4
H_238U = 9.37e-5;
H_235U = 5.69e-4;
H_232Th = 2.69e-5;
H_40K = 2.79e-5;
%Radioisotope abundances, Schubert et al. (2001) Ch. 4
Uran = 21e-9; %[U] = 21 ppb
C_238U = 0.9927 * Uran;
C_235U = 0.0072 * Uran;
C_40K = 1.28 * Uran;
C_232Th = 4.01 * Uran;
%Decay constants (1/yr), from McNamara and van Keken (2000)
l_238U = 0.155e-9;
l_235U = 0.985e-9;
l_232Th = 0.0495e-9;
l_40K = 0.555e-9;

%heat production (J/m^3/yr)
Q = (C_238U * H_238U * exp(l_238U*(4.6e9 -t)) + C_235U*H_235U*exp(l_235U*(4.6e9 -t)) + ...
    C_232Th * H_232Th * exp(l_232Th * (4.6e9 -t)) + C_40K * H_40K * exp(l_40K * (4.6e9 -t)))*...
    rho_m*3.15569e7;

[qm,Db,uc,Ra,num]= heatflux(T_water(1),f_water,Rp,Rc,g,avgfact,rho_m,avgP);

%spreading rate (m^2/s)
S =  2 * Lridge * uc; 

[rsub,~] = regassing(qm,f_water,Rp,m0,Mm);
[~,rmor,~,~] = degassing(T_water(1),Db,qm,f_water,Rp,g,avgfact);

%temperature derivative
deriv(1) = (-3 * Rp^2 * qm * 3.15569e7 + Q * (Rp^3 - Rc^3)) / (rho_m * cp * (Rp^3 - Rc^3)); %temperature per year
%water derivative
deriv(2) = (- rmor + rsub) * S * 3.15569e7; %mass fraction of water per year
% waterest = deriv(2) * 1e7 + T_water(2);
% if waterest > m0
%     deriv(2) = deriv2est;
% end



end
