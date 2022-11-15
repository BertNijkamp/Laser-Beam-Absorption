function [Pout,kout,eout,absorbed,R,Pinters] = Refraction3D(P,k,Pp,d,RI0,RI1,e,medium,stepsize,absorp)
% This function calculates the one or two new laser position(s) based on reflection and refraction,
% as well as their directions and relative energies and the total amount of absorbed energy.
% When entering or leaving more than one particle, the new surface normal is calculated based
% on how close the new position is to the particles.


Np = size(Pp,1);  % Number of particles the new ray position is inside of
r = d/2;

% Determine refractive index ratio, use medium 1 as reverser
if medium == 0
    ratio = RI0/RI1;
    reverser = 1;
else
    ratio = RI1/RI0;
    reverser = -1;
end

% Matrix used to calculate new coordinates
A = [k;P-Pp];
if Np > 1
    D = arrayfun(@(x)r(x)^2-sum(cellfun(@(n)det(A([1 x+1],n))^2,{1:2,[1 3],2:3})),1:Np)';
    Pinters = cell2mat(arrayfun(@(x)P-k*(sum(prod(A([1 x+1],:)))+reverser*sqrt(D(x))),(1:Np)','UniformOutput',false));
    [~,IntersParticle] = min(arrayfun(@(x)norm(P-Pinters(x,:)),1:Np));
    Pinters = Pinters(IntersParticle,:);
else
    D = r^2-sum(cellfun(@(n)det(A(:,n))^2,{1:2,[1 3],2:3}));
    Pinters = P-k*(sum(prod(A))+reverser*sqrt(D));
end

% Portion of every particle for resulting surface normal
if Np > 1
    portion = arrayfun(@(x)(1-(norm(Pinters-Pp(x,:))-r(x))/(.5*stepsize))/Np,1:Np)';
    portion(portion<0) = 0;
    portion(IntersParticle) = 1-sum(portion)+portion(IntersParticle);
else
    portion = 1;
end

normals = reverser*(Pinters-Pp)./arrayfun(@(x)norm(Pinters-Pp(x,:)),1:Np)';
normal = sum(normals.*portion,1);
normal = normal/norm(normal);

cosA = dot(-k,normal);
kr = k+normal*2*cosA;
travelled = norm(Pinters-P);

% For total internal reflection in particle
if medium==1 && abs(acos(cosA))>=asin(1/ratio)
    kout = kr;
    absorbed = (1-exp(-absorp*stepsize*10^6))*e;
    eout = e-absorbed;
    R = 1;

% For no total internal reflection
else
    cosB = sqrt(1-ratio^2*(1-cosA^2));
    kt = k*ratio+normal*(ratio*cosA-cosB);
    kout = [kt;kr];
    R = (((ratio*cosA-cosB)/(ratio*cosA+cosB))^2+((ratio*cosB-cosA)/(ratio*cosB+cosA))^2)/2;
    if medium == 0
        absorbed = (1-exp(-absorp*(stepsize-travelled)*10^6))*e*(1-R);
        eout = [e*(1-R)-absorbed;e*R];
    else
        absorbed1 = (1-exp(-absorp*travelled*10^6))*e;
        absorbed2 = (1-exp(-absorp*(stepsize-travelled)*10^6))*(e-absorbed1)*R;
        absorbed = absorbed1+absorbed2;
        eout = [(e-absorbed1)*(1-R);(e-absorbed1)*R-absorbed2];
    end
end
Pout = Pinters+kout*(stepsize-travelled);

end