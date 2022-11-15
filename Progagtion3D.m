clear,clc

SymmetricWalls = 0;
ReflectiveBottom = 0;
ShowLaser = 0;
ShowAbsorption = 0;
OnlyShowAbsorption = 0;

NumberOfRays = 10000;

% Load ParticleList
ParticleList = table2array(readtable('PG_60_1.txt'));
ParticleList(:,4) = ParticleList(:,4)*2;
ParticleList = [ParticleList';zeros(1,size(ParticleList,1))]';

% Determine bed dimensions
MinX = min(ParticleList(:,1)-ParticleList(:,4)./2);
MaxX = max(ParticleList(:,1)+ParticleList(:,4)./2);
MinY = min(ParticleList(:,2)-ParticleList(:,4)./2);
MaxY = max(ParticleList(:,2)+ParticleList(:,4)./2);
MinZ = min(ParticleList(:,3)-ParticleList(:,4)./2);
MaxZ = max(ParticleList(:,3)+ParticleList(:,4)./2);
BedX = MaxX - MinX;
BedY = MaxY - MinY;

LaserStepSize = min(ParticleList(:,4))/10;
P_Absorption = 0.0274;              % 1/µm for 532 nm wavelength (PS)   (Hejmady)
P_RefractiveIndex = 1.5997;         % Particle (medium 1) (PS)          (Sultanova)
A_RefractiveIndex = 1;              % Air (medium 0)
AbsorbedThr = 0.01/NumberOfRays;    % Stop propagating when energy of ray is 1%

% Storage decimals for multiple particle numbers
decimals = ceil(log10(size(ParticleList,1)));

% Start rays divided over circle
SpotRadius = 200e-6;  %m
LaserAngle = 0;       %degrees
StartingDirection = [0 tand(LaserAngle) -1];
Centre = [(MaxX+MinX)/2+((MaxZ-MinZ)*1.15-(MaxZ-MinZ)/2)*StartingDirection(1)/StartingDirection(3) ...
          (MaxY+MinY)/2+((MaxZ-MinZ)*1.15-(MaxZ-MinZ)/2)*StartingDirection(2)/StartingDirection(3)];

% Coords of laser beam
CircularCoords = PointsInCircle(Centre(1),Centre(2),SpotRadius,NumberOfRays);

% Energy of every ray
Energies = GaussianIntensity(CircularCoords,Centre,SpotRadius);

% Scale ellipse to circle in case of laser at angle
CircularCoords = (CircularCoords-Centre).*[norm(StartingDirection([1 3]))/abs(StartingDirection(3)) ...
                 norm(StartingDirection([2 3]))/abs(StartingDirection(3))]+Centre;

% LaserList variable numbers
% [1 2 3 4  5  6  7   8        9               10          ]
% [X Y Z kx ky kz e medium particlenr dissipated/out-of-bed]

% Make laser list for every ray
LaserList = cell(1,NumberOfRays);
for RayCounter = 1:size(CircularCoords,1)
    LaserList{RayCounter} = [CircularCoords(RayCounter,1) CircularCoords(RayCounter,2) (MaxZ-MinZ)*1.15+MinZ ...
                             StartingDirection Energies(RayCounter) 0 0 0];
    % Make direction a unit vector
    LaserList{RayCounter}(1,4:6) = LaserList{RayCounter}(1,4:6)/norm(LaserList{RayCounter}(1,4:6));
end

% Update bed size based on laser starting position
MinX = min([MinX Centre(1)-SpotRadius]);
MaxX = max([MaxX Centre(1)+SpotRadius]);
MinY = min([MinY Centre(2)-SpotRadius]);
MaxY = max([MaxY Centre(2)+SpotRadius]);

% Keep calculating positions untill all rays are dissipated/left the bed
while ~all(cellfun(@(x)x(end,10),LaserList))
    for RayNr = 1:size(LaserList,2)
        NewStep = size(LaserList{RayNr},1)+1;

        % Skip if ray is dissipated or out of bed
        if LaserList{RayNr}(NewStep-1,10) == 1
            continue
        end

        % Calculate next step of current ray
        LaserList{RayNr}(NewStep,1:3) = LaserList{RayNr}(NewStep-1,1:3)+LaserStepSize*LaserList{RayNr}(NewStep-1,4:6);

        % Test for change in medium
        [LaserList{RayNr}(NewStep,8),InParticleList] = InParticle3D(LaserList{RayNr}(NewStep,1:3),ParticleList);
        if size(InParticleList,2) > 1
            LaserList{RayNr}(NewStep,9) = MatrixToNumber(InParticleList,decimals);
        else
            LaserList{RayNr}(NewStep,9) = InParticleList;
        end
        if LaserList{RayNr}(NewStep,8) ~= LaserList{RayNr}(NewStep-1,8)

        % Steps if change in medium
            if LaserList{RayNr}(NewStep,8) == 1
                % If going into particle, use particlenr of new step
                ParticleNr = InParticleList;
            else
                % If going out of particle, use particlenr of old step
                ParticleNr = NumberToMatrix(LaserList{RayNr}(NewStep-1,9),decimals);
            end
            PartNum = size(ParticleNr,2);
            [Pout,kout,eout,absorbed] = Refraction3D(LaserList{RayNr}(NewStep-1,1:3),LaserList{RayNr}(NewStep-1,4:6), ...
                ParticleList(ParticleNr,1:3),ParticleList(ParticleNr,4),A_RefractiveIndex,P_RefractiveIndex, ...
                LaserList{RayNr}(NewStep-1,7),LaserList{RayNr}(NewStep-1,8),LaserStepSize,P_Absorption);
            if ~isreal(Pout)
                error('Imaginary results')
            end
            
            % Add absorbed energy to particle
            ParticleList(ParticleNr,5) = ParticleList(ParticleNr,5) + absorbed/PartNum;
            
            % Transmitted (or total internal reflected) rays
            LaserList{RayNr}(NewStep,1:3) = Pout(1,:);
            LaserList{RayNr}(NewStep,4:6) = kout(1,:);
            LaserList{RayNr}(NewStep,7) = eout(1);
            if size(Pout,1) == 1
                LaserList{RayNr}(NewStep,8) = 1;
                LaserList{RayNr}(NewStep,9) = MatrixToNumber(ParticleNr,decimals);
            end
            if LaserList{RayNr}(NewStep,7) < AbsorbedThr
                % If e is below threshold, set dissipated to 1 and add leftover energy to particle
                LaserList{RayNr}(NewStep,10) = 1;
                ParticleList(ParticleNr,5) = ParticleList(ParticleNr,5)+LaserList{RayNr}(NewStep,7)/PartNum;
            end
            
            % Reflected rays (skip for total internal reflection)
            if size(Pout,1) == 2
                % Create new ray
                LaserList{end+1}(1,1:3) = Pout(2,:); %#ok<SAGROW>
                LaserList{end}(1,4:6) = kout(2,:);
                LaserList{end}(1,7) = eout(2);
                [LaserList{end}(1,8),InParticleList] = InParticle3D(Pout(2,:),ParticleList);
                LaserList{end}(1,9) = MatrixToNumber(InParticleList,decimals);
                if LaserList{end}(1,7) < AbsorbedThr
                    % If e is below threshold, set dissipated to 1 and add leftover energy to particle
                    LaserList{end}(1,10) = 1;
                    ParticleList(ParticleNr,5) = ParticleList(ParticleNr,5)+LaserList{end}(1,7)/PartNum;
                else
                    LaserList{end}(1,10) = 0;
                end
            end
        else

        % Steps if no change in medium
            LaserList{RayNr}(NewStep,4:6) = LaserList{RayNr}(NewStep-1,4:6);
            if LaserList{RayNr}(NewStep,8) == 1
                ParticleNr = InParticleList;
                PartNum = size(ParticleNr,2);
                absorbed = (1-exp(-P_Absorption*LaserStepSize*10^6))*LaserList{RayNr}(NewStep-1,7);
                ParticleList(ParticleNr,5) = ParticleList(ParticleNr,5) + absorbed/PartNum;
                LaserList{RayNr}(NewStep,7) = LaserList{RayNr}(NewStep-1,7) - absorbed;
                if LaserList{RayNr}(NewStep,7) < AbsorbedThr
                    % If e is below threshold, set dissipated to 1 and add leftover energy to particle
                    LaserList{RayNr}(NewStep,10) = 1;
                    ParticleList(ParticleNr,5) = ParticleList(ParticleNr,5)+LaserList{RayNr}(NewStep,7)/PartNum;
                end
            else
                LaserList{RayNr}(NewStep,7) = LaserList{RayNr}(NewStep-1,7);
            end
        end
        
        % Check if outside of vertical borders, teleport if symmetricwalls is on
        Teleported = 0;
        Reflected = 0;
        if LaserList{RayNr}(NewStep,1) < MinX
            if SymmetricWalls == 1
                LaserList{RayNr}(NewStep,1) = LaserList{RayNr}(NewStep,1) + BedX;
                Teleported = 1;
            else
                LaserList{RayNr}(NewStep,10) = 1;
            end
        elseif LaserList{RayNr}(NewStep,1) > MaxX
            if SymmetricWalls == 1
                LaserList{RayNr}(NewStep,1) = LaserList{RayNr}(NewStep,1) - BedX;
                Teleported = 1;
            else
                LaserList{RayNr}(NewStep,10) = 1;
            end
        end
        if LaserList{RayNr}(NewStep,2) < MinY
            if SymmetricWalls == 1
                LaserList{RayNr}(NewStep,2) = LaserList{RayNr}(NewStep,2) + BedY;
                Teleported = 1;
            else
                LaserList{RayNr}(NewStep,10) = 1;
            end
        elseif LaserList{RayNr}(NewStep,2) > MaxY
            if SymmetricWalls == 1
                LaserList{RayNr}(NewStep,2) = LaserList{RayNr}(NewStep,1) - BedY;
                Teleported = 1;
            else
                LaserList{RayNr}(NewStep,10) = 1;
            end
        end
        
        % Check if outside of Z borders
        if LaserList{RayNr}(NewStep,3) > (MaxZ-MinZ)*1.15
            LaserList{RayNr}(NewStep,10) = 1;
        elseif LaserList{RayNr}(NewStep,3) < MinZ
            % Reflect if reflective bottom is on
            if ReflectiveBottom == 1
                LaserList{RayNr}(NewStep,3) = -LaserList{RayNr}(NewStep,3);
                LaserList{RayNr}(NewStep,6) = -LaserList{RayNr}(NewStep,6);
                Reflected = 1;
            else
                LaserList{RayNr}(NewStep,10) = 1;
            end
        end
        
        % Recheck whether in particle if teleported/reflected
        if Teleported==1 || Reflected==1
            [LaserList{RayNr}(NewStep,8),InParticleList] = InParticle3D(LaserList{RayNr}(NewStep,1:3),ParticleList);
            LaserList{RayNr}(NewStep,9) = MatrixToNumber(InParticleList,decimals);
            Teleported = 0;
            Reflected = 0;
        end
    end

    clc,disp(['Calculating laser paths, ' num2str(sum(~cellfun(@(x)x(end,10),LaserList))) ' rays left...'])
end

% Draw particles and laser
if ShowLaser == 1
    clc,disp('Drawing particles and lasers...')
    DrawParticles3D(ParticleList)
    hold on
    xlim([MinX MaxX])
    ylim([MinY MaxY])
    zlim([MinZ LaserList{1}(1,3)])
    for RayNr = 1:size(LaserList,2)
        for Step = 1:size(LaserList{RayNr},1)-1
            % Skip if line teleports to other side
            if SymmetricWalls==1 && norm(LaserList{RayNr}(Step+1,1:3)-LaserList{RayNr}(Step,1:3)) > 2*LaserStepSize
                continue
            end
            MaxColour = max(Energies);
            line(LaserList{RayNr}(Step:Step+1,1),LaserList{RayNr}(Step:Step+1,2),LaserList{RayNr}(Step:Step+1,3),...
                'Color',[1 0 0 .75*LaserList{RayNr}(Step,7)/MaxColour],'LineWidth',2)
        end
    end
    view(250,40)
end

% Draw particles with absorption
if ShowAbsorption == 1
    clc,disp('Drawing particles with absorption...')
    if OnlyShowAbsorption == 1
        DrawAbsorption3D(ParticleList(ParticleList(:,5)~=0,:))
        xlim([MinX MaxX])
        ylim([MinY MaxY])
        zlim([MinZ MaxZ])
    else
        DrawAbsorption3D(ParticleList)
    end
    set(gca,'color','none'); 
end

% Result
disp([num2str((1-sum(ParticleList(:,5)))*100) '% of the energy left the bed.'])
ParticleTable = array2table(ParticleList);
writetable(ParticleTable,'Absorption_60_1.txt','WriteVariableNames',false)
