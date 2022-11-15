function [] = DrawAbsorption3D(ParticleList)
% This function renders the spherical particles with a colour based on their absorption.

islog = 0;      % 0=linear 1=logarithmic

% Calculate boundaries
MinX = min(ParticleList(:,1)-ParticleList(:,4)./2);
MaxX = max(ParticleList(:,1)+ParticleList(:,4)./2);
MinY = min(ParticleList(:,2)-ParticleList(:,4)./2);
MaxY = max(ParticleList(:,2)+ParticleList(:,4)./2);
MinZ = min(ParticleList(:,3)-ParticleList(:,4)./2);
MaxZ = max(ParticleList(:,3)+ParticleList(:,4)./2);
if ~islog
    Absorption = ParticleList(:,5);
else
    Absorption = log(ParticleList(:,5));
    Absorption(isinf(Absorption)) = min(Absorption(~isinf(Absorption)));
end
MinAbsorption = min(Absorption);
MaxAbsorption = max(Absorption);
colour = 1-(Absorption-MinAbsorption)/(MaxAbsorption-MinAbsorption);
colour(colour>1) = 1;

% Draw particles
clf
plot3(0,0,0)
hold on
axis equal
axis([MinX MaxX MinY MaxY MinZ MaxZ])

for count = 1:size(ParticleList,1)
    X = ParticleList(count,1);
    Y = ParticleList(count,2);
    Z = ParticleList(count,3);
    R = ParticleList(count,4)/2;
    [x,y,z] = sphere;
    surf(R*x+X,R*y+Y,R*z+Z,'LineStyle','none','FaceAlpha',1,'FaceColor',[1 colour(count)*ones(1,2)])
    xlabel('x [m]'),ylabel('y [m]'),zlabel('z [m]')
end

light('Position',[-1 -.5 2],'Color',[.8 .8 .8])
light('Position',[1 1 -1],'Color',[.2 .2 .2])
material dull

end