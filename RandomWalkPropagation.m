clear,clc

% Bed properties
BedX = 1000e-6; %m
BedY = 1000e-6; %m
BedZ = 500e-6; %m

% Laser properties
NumberOfRays = 10000;
AbsorpP = 27400;                          %1/m
RefrIndexP = 1.5997;
RefrIndexA = 1;
AbsorbedThr = 0.001/NumberOfRays;
ShowLaser = 1;
LaserType = 1;                            % 0=solid colour,     1=lin scale, 2=log scale
GridType = 0;                             % 0=random positions, 1=grid
LaserAngle = 0;                           %degrees
LaserList = cell(1,NumberOfRays);

% Particle properties
AvDiameter = 60e-6;                       %m
sigma = 0;                                % Constant Diameter
MinDiameter = 30e-6;                      %m
MaxDiameter = 90e-6;                      %m
VolumeFraction = 0.6;
ParticleList = zeros(50*NumberOfRays,5);  % Allocate enough space for particle list
LastRow = 1;

% Starting values at top of bed
for RayNr = 1:NumberOfRays
    if GridType
        p = [BedX/sqrt(NumberOfRays)*(rem(RayNr-1,sqrt(NumberOfRays))+.5) ... 
             BedY/sqrt(NumberOfRays)*(ceil(RayNr/sqrt(NumberOfRays))-.5) BedZ];
    else
        p = [rand*BedX rand*BedY BedZ];
    end
    k = [0 tand(LaserAngle) -1];
    e = 1/NumberOfRays;
    medium = 0;
    pp = [NaN NaN NaN];
    D = NaN;
    L = FirstAirLength(AvDiameter);
    LaserList{RayNr} = [p k e medium pp D L 0];
end

% Simulate random walk propagation for each ray
while ~all(cellfun(@(x)x(end,14),LaserList))
    for RayNr = 1:size(LaserList,2)
        NewStep = size(LaserList{RayNr},1)+1;
        
        % Skip if ray is dissipated or out of bed
        if LaserList{RayNr}(NewStep-1,14) == 1
            continue
        end
        
        p_old = LaserList{RayNr}(NewStep-1,1:3);
        k_old = LaserList{RayNr}(NewStep-1,4:6);
        e_old = LaserList{RayNr}(NewStep-1,7);
        L_old = LaserList{RayNr}(NewStep-1,13);
        p_new = p_old+k_old*L_old;
        
        % Going into particle:
        if LaserList{RayNr}(NewStep-1,8) == 0
            
            % Test for going out of bed
            if p_new(3)>BedZ
                LaserList{RayNr}(NewStep,:) = [p_new NaN*ones(1,10) 1];
                continue
            end
            
            % Position of new particle
            BelowMaxZ = 0;
            counter = 0;
            while ~BelowMaxZ && counter<100
                D_new = Diameter(AvDiameter,sigma,MinDiameter,MaxDiameter);     % random value
                r = D_new/2*sqrt(rand);                                         % random value
                theta = rand*2*pi;                                              % random value
                h = [0 -k_old(3) k_old(2)];
                pp_new = p_new+k_old*sqrt(D_new^2/4-r^2)-r*(h*cos(theta)+cross(k_old,h)*sin(theta));
                if pp_new(3)<BedZ-D_new/2
                    BelowMaxZ = 1;
                end
                counter = counter+1;
            end
            
            % If failed to place particle, treat as going out of bed
            if counter == 100
                LaserList{RayNr}(NewStep,:) = [p_new NaN*ones(1,10) 1];
                continue
            end
            
            % Test for dissipated
            if e_old<AbsorbedThr
                LaserList{RayNr}(NewStep-1,14) = 1;
                LaserList{RayNr}(NewStep,:) = [p_new NaN*ones(1,10) 1];
                ParticleList(LastRow,:) = [pp_new D_new e_old];
                LastRow = LastRow+1;
                continue
            end
            
            eta = RefrIndexA/RefrIndexP;
            n = p_new-pp_new;
            n = n/norm(n);
            cosA = dot(-k_old,n);
            cosB = sqrt(1-eta^2*(1-cosA^2));
            R = (((eta*cosA-cosB)/(eta*cosA+cosB))^2+((cosA-eta*cosB)/(cosA+eta*cosB))^2)/2;
            
            % Reflected ray
            k_newR = k_old+2*n*cosA;
            e_newR = e_old*R;
            L_newR = AirLength(VolumeFraction,AvDiameter);
            LaserList{end+1} = [p_new k_newR e_newR 0 [NaN NaN NaN] NaN L_newR 0]; %#ok<SAGROW>
            
            % Transmitted ray
            k_new = eta*k_old+n*(eta*cosA-cosB);
            e_new = e_old-e_newR;
            L_new = D_new*cosB;
            LaserList{RayNr}(NewStep,:) = [p_new k_new e_new 1 pp_new D_new L_new 0];
            if e_new<AbsorbedThr
                LaserList{RayNr}(NewStep,14) = 1;
                ParticleList(LastRow,:) = [pp_new D_new e_new];
                LastRow = LastRow+1;
            end
            
        % Going out of particle:
        else
            pp_old = LaserList{RayNr}(NewStep-1,9:11);
            D_old = LaserList{RayNr}(NewStep-1,12);
            eta = RefrIndexP/RefrIndexA;
            n = pp_old-p_new;
            n = n/norm(n);
            cosA = dot(-k_old,n);
            e_abs = e_old*exp(-AbsorpP*L_old);

            % For total internal reflection:
            if abs(acos(cosA)) >= asin(1/eta)
                LaserList{RayNr}(NewStep,:) = [NaN*ones(1,13) 1];
                ParticleList(LastRow,:) = [pp_old D_old e_old];
                LastRow = LastRow+1;

            % For no total internal reflection:
            else
                ParticleList(LastRow,:) = [pp_old D_old e_old-e_abs];
                LastRow = LastRow+1;
                cosB = sqrt(1-eta^2*(1-cosA^2));
                R = (((eta*cosA-cosB)/(eta*cosA+cosB))^2+((cosA-eta*cosB)/(cosA+eta*cosB))^2)/2;

                % Reflected ray
                k_newR = k_old+2*n*cosA;
                e_newR = e_abs*R;
                if e_newR>AbsorbedThr
                    LaserList{end+1} = [p_new k_newR e_newR 1 pp_old D_old L_old 0]; %#ok<SAGROW>
                else
                    ParticleList(LastRow,:) = [pp_old D_old e_newR];
                    LastRow = LastRow+1;
                end

                % Transmitted ray
                k_new = eta*k_old+n*(eta*cosA-cosB);
                e_new = e_abs-e_newR;
                L_new = AirLength(VolumeFraction,AvDiameter);
                LaserList{RayNr}(NewStep,:) = [p_new k_new e_new 0 [NaN NaN NaN] NaN L_new 0];
                if e_new<AbsorbedThr
                    LaserList{RayNr}(NewStep,14) = 1;
                end
            end
        end
        
    end
end

ParticleList(LastRow:end,:) = [];   % Remove leftover empty spots
NumberOfLayers = round((BedZ-AvDiameter)/(AvDiameter*sqrt(6)/3)+1);
LayerHeight = (BedZ-AvDiameter)/(NumberOfLayers-1);
EnergyPerLayer = zeros(NumberOfLayers,1);
for counter = 1:size(ParticleList,1)
    LayerNr = ceil((BedZ-AvDiameter/2-ParticleList(counter,3))/LayerHeight);
    if ParticleList(counter,1) < 0 || ParticleList(counter,1) > BedX
        continue
    elseif ParticleList(counter,2) < 0 || ParticleList(counter,2) > BedY
        continue
    elseif ParticleList(counter,3) < 0
        continue
    end
    EnergyPerLayer(LayerNr) = EnergyPerLayer(LayerNr)+ParticleList(counter,5);
end

% Plot absorption
semilogx(EnergyPerLayer,1:length(EnergyPerLayer))
set(gca,'YDir','reverse')
xlabel('Total E_{abs} in layer')
ylabel('Layer nr [-]')
title('Energy absorption in semi-random walk')

% Render laser rays
if ShowLaser == 1
    figure; view(3)
    hold on
    axis equal
    axis([min(cellfun(@(x) min(x(:,1)),LaserList)) max(cellfun(@(x) max(x(:,1)),LaserList)) ...
          min(cellfun(@(x) min(x(:,2)),LaserList)) max(cellfun(@(x) max(x(:,2)),LaserList)) ...
          min(cellfun(@(x) min(x(:,3)),LaserList)) max(cellfun(@(x) max(x(:,3)),LaserList))])
    for RayNr = 1:size(LaserList,2)
        for Step = 1:size(LaserList{RayNr},1)-1
            if LaserType == 0 % Solid colour
                line(LaserList{RayNr}(Step:Step+1,1),LaserList{RayNr}(Step:Step+1,2),LaserList{RayNr}(Step:Step+1,3),...
                    'Color',[1 0 0],'LineWidth',2)
            elseif LaserType == 1 % Linear scale
                line(LaserList{RayNr}(Step:Step+1,1),LaserList{RayNr}(Step:Step+1,2),LaserList{RayNr}(Step:Step+1,3),...
                    'Color',[1 0 0 LaserList{RayNr}(Step,7)*NumberOfRays],'LineWidth',2)
            elseif LaserType == 2 % Logarithmic scale
                line(LaserList{RayNr}(Step:Step+1,1),LaserList{RayNr}(Step:Step+1,2),LaserList{RayNr}(Step:Step+1,3),...
                    'Color',[1 0 0 max(1-log(LaserList{RayNr}(Step,7))/log(AbsorbedThr),0)],'LineWidth',1)
            end
        end
    end
    view(20,20)
end

% Calculation of first distance in air
function [l] = FirstAirLength(D)
    mu = -1.48462+log(D);
    sigma = 0.919366;
    l = lognrnd(mu,sigma);
end

% Calculation of distance in air
function [l] = AirLength(VF,D)
    mu = -2.7136*VF+0.4278+log(D);
    sigma = -0.5846*VF+1.212;
    l = lognrnd(mu,sigma);
end

% Calculation of particle size
function [D] = Diameter(mu,sigma,min,max)
    D = rand*sigma+mu;
    while D<min || D>max
        D = rand*sigma+mu;
    end
end
