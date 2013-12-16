function [XV] =SetUpInitDistro(vF,vSpread)

%This detuning is in angular frequency.
%Constants
kB = 1.38e-23;
mass = 105*1.6605e-27;
T = 3.2;
% vF = 78;
CellToMotDist = 0.4;
c = 3e8;
ApertureSize = .005;

%Initialize molecules.
% YAGFire = 3e-3;
YAGFire = 0;
CellEmptyingTime= 2e-3;
 NumberOfMolecules = 1e6;
% NumberOfMolecules = 10000;

XV = zeros(10,NumberOfMolecules); %[x y z vx vy vz CellExit vzInit State NumberOfTrappableMoleculesIndex]
XV(9,:) = ones(1,NumberOfMolecules);

%fGaussian = (mass/(2 *pi*kB*T))^(1/2)*exp(-mass*(v - vf).^2/(2*kB*T));
%Velocity Distribution

% XV(6,:) = normrnd(0,sqrt(kB*T/mass),1,NumberOfMolecules)+vF;

XV(6,:) = normrnd(0,vSpread,1,NumberOfMolecules)+vF;

XV(8,:) = XV(6,:);%8th row is the original velocities.


XV(4,:) = normrnd(0,sqrt(kB*T/mass),1,NumberOfMolecules); %generate thermal velocities for x and y
XV(5,:) = normrnd(0,sqrt(kB*T/mass),1,NumberOfMolecules);
XV = XV(:,XV(6,:)>0); %only allow the ones that are going forward
NumberOfMolecules = size(XV,2);
XV(7,:) = -CellEmptyingTime*log(1-rand(1,NumberOfMolecules)) + YAGFire; %Cell Exit time.
% XV(7,:) = XV(7,:) + .01255 +.0075;
% XV(7,:) = XV(7,:) -.00006582.*XV(6,:);
% XV(7,:) = XV(7,:) -.00012.*XV(6,:);
% XV(7,:) = XV(7,:)- CellToMotDist./(XV(6,:)); %Cell Exit time.

XV = XV(:,XV(7,:)>0); %only allow positive exit times
NumberOfMolecules = size(XV,2);
%Assume Point Source.
radius=rand(1,NumberOfMolecules);
theta = 2*pi*rand(1,NumberOfMolecules);
XV(3,:) = -XV(6,:).*XV(7,:);%zinit = -vzint*texit i.e. that at the cell exit time the molecule is at z = 0.
XV(2,:) = -XV(5,:).*XV(7,:) + ApertureSize*cos(theta).*sqrt(radius);
XV(1,:) = -XV(4,:).*XV(7,:)+ ApertureSize*sin(theta).*sqrt(radius);

% XV = XV(:,abs(XV(4,:)./XV(6,:))<=0.1&abs(XV(5,:)./XV(6,:))<=0.1);%throw
% out ones that are fast in x or y compared to z


save('XVOrigTemp','XV')
end
