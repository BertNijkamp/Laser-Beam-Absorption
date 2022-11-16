function [I] = GaussianIntensity(Coords,centre,SpotRadius)
% This function calculates the intensity for all rays in a circular pattern. The total intensity is 1.

% Distance from every point to the centre
d = arrayfun(@(n)norm([Coords(n,1)-centre(1) Coords(n,2)-centre(2)]),(1:size(Coords,1))');

% Intensity per ray
I = exp(-2*d.^2/SpotRadius^2)*2/(size(Coords,1)*(1-exp(-2)));
I = I/sum(I);

end
