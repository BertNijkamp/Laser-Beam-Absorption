function [diamOut] = DiameterIncrease(diam,D)
% This function calculates the increase in particle diameter due to two particles overlapping.
% The particles need to have equal diameters.

if D < diam
    Drel = D/diam;
    A = (Drel^3+64+8*sqrt(2*Drel^3+64))^(1/3);
    factor = (Drel^2/A-Drel+A)/4;
    diamOut = diam*factor;
else
    diamOut = diam;
end

end