function [diamOut] = DiameterIncrease(diam,D)
% This function calculates the increase in particle diameter due to two particles overlapping
% Warning: particles need to have equal diameters

if D < diam
    Drel = D/diam;
    A = (Drel^3+8+4*sqrt(Drel^3+4))^(1/3);
    factor = (Drel^2/A-Drel+A)/2;
    diamOut = diam*factor;
else
    diamOut = diam;
end

end