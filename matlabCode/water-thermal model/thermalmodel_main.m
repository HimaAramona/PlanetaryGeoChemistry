% This script solves the model of Sandu et al. (2011) for the mantle water
% content and thermal history. 

% I'm using the ODE solver to solve the mantle heat balance equation and
% the water exchange between surface and mantle. 
clear;clc;

%Earth parameters
MEarth = 5.97e24; %kg
REarth = 6371e3; %km
%ocean mass (kg)
ocean_kg = 1.39e21;
%Gravitational acceleration
G = 6.67384e-11; %m^3/kg/s^2

Mp = input('Planet Mass (in Earth Masses): ')*MEarth; % kg
Rp = input('Planet Radius (in Earth Radii): ')*REarth;%km
Rc = input('Core Radius (in Earth Radii): ')*REarth;%km
Xcore = .3259;
Xmantle = 1-.3259;
XH2O = input('Water mass fraction: ');
filename = input('Filename to save variables (ending with .mat): ');

Mm = Xmantle * Mp;
Mcore = Xcore * Mp;
MH2O = XH2O * Mp;
g = G * Mp / Rp^2;
alpha = 2e-5;
cp = 1.2e3;
%mantle density (kg/m3)
rho_m = 3 * Mm / (4 * pi() * (Rp^3 - Rc^3));


%need normalizing factor to find potential temperature from average mantle
%temperature for Db/q calculations
r = (Rc:1e3:Rp)';
top = trapz(r,r.^2 .* (1+alpha * g * (Rp - r)/cp));
bottom = trapz(r,r.^2);
avgfact = top/bottom;
topP = trapz(r,r.^2 .* (rho_m * g * (Rp - r)));
bottomP = trapz(r,r.^2);
avgP = topP/bottomP;

%initial conditions
T0 = 2520*avgfact; %initial temperature in K
m0 = MH2O; %intial water abundance in kg

%solve ODE for mantle temperature and water content evolution
%t is time in years, integration goes to 10 Gyr
%T_water(:,1) = upper mantle homologous temperature
%T_water(:,2) = mass of water in mantle (kg)
options = odeset('NonNegative',1:2,'RelTol',1e-4);
[t,T_water] = ode45(@(t,T_water) derivatives(t,T_water,Mm,Rp,Rc,g,avgfact,...
    m0,avgP),[0 10e9],[T0 m0],options);

% mid-ocean ridge length = 1.5x planet circumference (m)
Lridge = 2 * pi() * Rp * 1.5 ;
%Mass of the mantle (kg)
% Mmantle = 4.06e24; 
f_water = T_water(:,2)/ Mm;


qm = zeros(length(t),1);Db = zeros(length(t),1);uc = zeros(length(t),1);
rsub = zeros(length(t),1); Dhyd = zeros(length(t),1);
rmor = zeros(length(t),1);Dmelt = zeros(length(t),1);avgfmelt = zeros(length(t),1);
avgXmelt = zeros(length(t),1);S = zeros(length(t),1);

%Then calculate all temperature dependent variables
for i = 1:length(t)
    [qm(i),Db(i),uc(i),Ra(i),num(i)] = heatflux(T_water(i,1),f_water(i),Rp,Rc,g,avgfact,rho_m,avgP);
    [rsub(i),Dhyd(i)] = regassing(qm(i),f_water(i),Rp,m0,Mm);
    [Dmelt(i),rmor(i),avgfmelt(i),avgXmelt(i)] = degassing(T_water(i,1),Db(i),qm(i),f_water(i),Rp,g,avgfact);
    S(i) = 2 * Lridge * uc(i);
end

save(filename,'t','T_water','f_water','qm','Db','uc','Dhyd','Dmelt','num','avgfmelt',...
    'avgXmelt','S','Rp','Rc','Mp','Xcore','Xmantle','XH2O','avgfact','Ra','rmor','rsub');

