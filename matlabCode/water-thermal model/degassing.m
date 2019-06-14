function [Dmelt,rmor,avgfmelt,avgXmelt] = degassing(temp,Db,qm,f_water,Rp,g,avgfact)

%This function calculates the degassing rate. To do this, we must first
%calculate the melt fraction in the melt zone, the weight fraction of water
%in the melt and the thickness of the melt zone. 

%bulk distribution of water between silicate and melt
D_H2O = 0.01;

%density
rho_m = 3.3e3;
beta = (3/2);
% parameters to calculate offset from solidus due to water
K = 43; %degrees C / wt%^gamma
gamma = (3/4);
alpha = 2e-5;
cp=1.2e3;
chi_d = 0.02;
km = 4.2;
% Rp = 6271e3;

%pressure (GPa) from surface to an arbitrary depth. 300 km is what Sandu et
%al. used. May need to modify this in future.
z = (0:1e3:300e3)'; %meters
y = zeros(length(z),1);
Xmelt = zeros(length(z),1);
Tprofile = zeros(length(z),1);

press = rho_m * g * z/1e9; %in GPa
for i = 1:length(z)
    if press(i) < 10e0;
        Tsolidus(i) = 1120.661 + 132.899 * press(i) - 5.104 * press(i).^2+273; %C
    elseif press(i) >= 10.0e0 && press(i) < 23.5e0;
        Tsolidus(i) = -1.092 * (press(i)-10)^2 + 32.39 * (press(i) - 10)+1939.251 + 273;
    elseif press(i) >= 23.5 e0;
        Tsolidus(i) = 26.53 * (press(i) - 23.5) + 2175. + 273;
    end
end

%liquidus fit from fig. 15 of Zhang & Herzberg 1994
Tliquidus = 1710 + 59.062 * press - 2.9121 * press.^2 + 0.072019 * press.^3 + 273;
Tsol_wet = Tsolidus - K * (100*f_water/D_H2O)^gamma;
Tliq_wet = Tliquidus - K * (100*f_water)^gamma;

%calculate temperature profile through the boundary layer and upper mantle
Tprofile(1) = 300;%K
%mantle potential temperature
Tp = temp / avgfact;
for i = 2:length(z)
    if z(i) < Db
        %boundary layer thermal conductive temperature profile
        Tprofile(i) = Tprofile(1) + z(i) * qm/km; %K
    else
%         mantle adiabat
        Tprofile(i) = Tp + Tp *(alpha*g*(z(i))/cp);
    end
end

y(1)=0;
Xmelt(1) = 0;
options = optimset('Display','off');

%try finding Fmelt using fzero or fsolve
for i = 2:length(z)
    
    if Tprofile(i) >= Tliq_wet(i)
        
        y(i) = 1.0;
                
    elseif Tprofile(i) < Tliq_wet(i) && Tprofile(i) >= Tsol_wet(i)
        
        h = @(f)abs(f-((Tprofile(i) - (Tsolidus(i) - K * (100*(f_water)/(D_H2O + f...
            *(1-D_H2O)))^gamma))/(Tliquidus(i) - Tsolidus(i)))^beta);
        y(i) = fsolve(h,0.1,options);
        
        if y(i) <= 0.0;
            y(i) = 0.0;
        elseif y(i) > 1.0;
            y(i) = 1.0;
        end
        
    else
        
        y(i) = 0;
        
    end
    
    if y(i) > 0e0;
        
        Xmelt(i) = (f_water/(D_H2O + y(i)*(1-D_H2O)));
        
    else
        
        Xmelt(i) = 0e0;
        
    end
    
end

%calculate the thickness of the melt layer
%the solidus/liquidus has an turning point at the mid-mantle, so it's
%possible to have a (partially) solid layer between two melt layers. The
%loop below finds the boundaries of the melt zones.
firstz = find(y,1,'first');
secondz = 1;
thirdz = 1;
fourthz = 1;
for i = firstz:length(z)
    if y(i) > 0.0e0;
        secondz = int16(i);
    end
    if i == length(z);
        break
    end
    if y(i+1) <= 0.0e0;
        break
    end
end
n = secondz + 1;
if sum(y(n:end)) > 0;
    y2 = y(n:end);
    thirdz = find(y2,1,'first')+n-1;
    fourthz = length(z);
end

% last = find(y,1,'last');

if isempty(firstz)
    firstz = int16(1);
end

Dmelt = z(secondz)-z(firstz) + z(fourthz) - z(thirdz);

r = Rp - z;

%find average water and melt fraction over the melt layer thickness
if Dmelt > 0e0;
    if fourthz > 1
        product1 = trapz(r(firstz:fourthz),y(firstz:fourthz).*r(firstz:fourthz).^2);
        product2 = trapz(r(firstz:fourthz),Xmelt(firstz:fourthz).*r(firstz:fourthz).^2);
        avgfmelt = 3* product1/(r(secondz)^3 - r(firstz)^3 + r(fourthz)^3 - r(thirdz)^3);
        avgXmelt = 3* product2/(r(secondz)^3 - r(firstz)^3 + r(fourthz)^3 - r(thirdz)^3);
    else
        product1 = trapz(r(firstz:secondz),y(firstz:secondz).*r(firstz:secondz).^2);
        product2 = trapz(r(firstz:secondz),Xmelt(firstz:secondz).*r(firstz:secondz).^2);
        avgfmelt = 3* product1/(r(secondz)^3 - r(firstz)^3);
        avgXmelt = 3* product2/(r(secondz)^3 - r(firstz)^3);
    end
    rmor = avgfmelt * avgXmelt * Dmelt * chi_d * rho_m;
else
    avgfmelt = 0e0;
    avgXmelt = 0e0;
    rmor = 0e0;
end


end





