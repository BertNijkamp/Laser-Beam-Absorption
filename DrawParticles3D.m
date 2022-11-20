function [] = DrawParticles3D(ParticleList)
% This function renders the spherical particles with a colour based on their diameter.

% Calculate boundaries
MinX = min(ParticleList(:,1)-ParticleList(:,4)./2);
MaxX = max(ParticleList(:,1)+ParticleList(:,4)./2);
MinY = min(ParticleList(:,2)-ParticleList(:,4)./2);
MaxY = max(ParticleList(:,2)+ParticleList(:,4)./2);
MinZ = min(ParticleList(:,3)-ParticleList(:,4)./2);
MaxZ = max(ParticleList(:,3)+ParticleList(:,4)./2);
MinD = min(ParticleList(:,4));
MaxD = max(ParticleList(:,4));

% No colours if radii are equal
if MinD == MaxD
    MakeColour = 0;
else
    MakeColour = 1;
end

% Draw particles
clf
view(axes,3)
hold on
axis equal
axis([MinX MaxX MinY MaxY MinZ MaxZ])

for count = 1:size(ParticleList,1)
    X = ParticleList(count,1);
    Y = ParticleList(count,2);
    Z = ParticleList(count,3);
    R = ParticleList(count,4)/2;
    if MakeColour == 1
    	colour = ((1-(2*R-MinD)/(MaxD-MinD))*ones(1,3))*0.7+0.2;
    else
        colour = [1 1 1];
    end
    [x,y,z] = sphere;
    surf(R*x+X,R*y+Y,R*z+Z,'LineStyle','none','FaceColor',colour)
end

xlabel('x [m]'),ylabel('y [m]'),zlabel('z [m]')
light('Position',[-1 -.5 2],'Color',[.8 .8 .8])
light('Position',[1 1 -1],'Color',[.2 .2 .2])
material dull

end
