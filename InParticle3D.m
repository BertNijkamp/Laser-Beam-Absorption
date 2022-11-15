function [in,nr] = InParticle3D(P,ParticleList)
% This function checks if the new laser position is inside any particle.
% It outputs yes/no and the number of the particle(s) it is in.

in = 0;
nr = [];

% Test whether current coordinate is inside any particle
for i = 1:size(ParticleList,1)
    distance = norm(P-ParticleList(i,1:3));
    if distance < ParticleList(i,4)/2
        in = 1;
        nr = [nr i]; %#ok<AGROW>
    end
end

% Set particle nr to 0 if not in particle
if in == 0
    nr = 0;
end

end