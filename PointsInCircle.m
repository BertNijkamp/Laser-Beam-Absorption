function [XY] = PointsInCircle(X,Y,r,N)
% This function randomly distributes N points in a circle with radius r and centre (X,Y).

rsquared = rand(N,1)*r.^2;
theta = rand(N,1)*2*pi;
XY = [sqrt(rsquared).*cos(theta)+X sqrt(rsquared).*sin(theta)+Y];

end