clear
%load('Sept25noacc.mat')
% load('simnoacc.mat'); %taken with 1.4*10^6 initial molecules
load('Aug_28_191329.mat'); %unslowed beam in mot region
motmcs = mcs;
motscope = scope;

load ('Nov_05_150456.mat');
cubemcs = mcs;
cubescope = scope;
cubeyag = .0012;
% FluorvTimeNoAcc = FluorvTime;

% angle = [0 -.02 -.04 -.06 -.08 -.1];


%  motDetuning = 1.6*detectVels*10^6; %Hz, compared to absolute frequency

gammaV0 = 2*pi*5*10^6;
kB = 1.38e-23;
mass = 105*1.6605e-27;
T = 3.2;
CellToMotDist = 0.40;
CellToCubeDist = .20;
CubeBeamSize = .002;
MOTsize = .010;
TrappingV = 10;
c = 3e8;
lambda1 = 648e-9;
lambda0 = 614e-9;
f0 = c/lambda0;
detectionPower = .005*10^-3; %for 1 cm beam
satIntensity0 = 2.3*10^-3; %per cm^2
s0 = detectionPower/satIntensity0;


v2angle = .01; %in rads, angle between v=1, 2 in slowing beam;
detectionAngle = .525; %in rads, angle from normal to beam;
% detectionAngle = 0;
% detectionAngle = 0; %in rads, angle from normal to beam;

% motDetuning = 0*10^6; %Hz, compared to absolute frequency
% detectVels = [50 60 70 80 90 100 110 120 130];
detectVels = [40 60];
%motDetuning = -1.6*detectVels*10^6;
%motDetuning = -1.6*detectVels*10^6*sin(detectionAngle); %Hz, compared to absolute frequency
%motDetuning = -f0*detectVels*sin(detectionAngle)/c./(1.+detectVels*sin(detectionAngle)/c);
motDetuning = -f0*detectVels*sin(detectionAngle)/c;
averagenum = [1 20];
% SweepTimeEnd = [8 9 10 11 12 13]*10^-3;
% vFast = [105 90 75 60 45];
% deltaV = [20 40 60];
vF = [80];
vSpread = [16];
for l = 1:size(detectVels,2) %go through velocity detunings
for k = averagenum(1):averagenum(2) % go through process multiple times for averaging purposes
    for i =1:size(vF,2) % go through different forward velocities
        for j = 1:size(vSpread,2) % go through different forward velocity spreads
            rng('shuffle');
            
            
            XV = SetUpInitDistro(vF(i),vSpread(j));
%             load('XVOrigTemp') %Gotten from SetupInitDistro.m
            %             NumberOfMolecules = 1.4e6;
            NumberOfMolecules = size(XV,2)
            % NumberOfMolecules = 100;
            
            
            
            
            %angle = .1; %angle between v=0/v=2 beams in radians, intersecting at MOT region
            v2BeamSize = .005; %radius
            
            hFine = 1e-6; %simulation time step
            hCoarse = 1e-4;
            SimTime = 100e-3;
            time = 0;
            vint = zeros(3,NumberOfMolecules);
            acc = zeros(3,NumberOfMolecules);
            SimIndex = 1;
            
            %tic
            %yint = HomeMade1DSplint(x,AccVsDetuning,Accpp);
            % FluorvTime = zeros(2,round(SimTime/hFine));
            NumberOfMolecules=size(XV,2);
            vStart = zeros(1,NumberOfMolecules);
            NumberOfTrappable = 0;
            LostToDelta = 0;
            
            %Frequency Sweep Parameters.

            vFast = 60;
            vSlow = 40;
            SweepTimeStart =5e-3;
            SweepTimeEnd = 9e-3;
            %             vSlow = vFast(j) - deltaV(k);
            
            SweepFrequencyStart = -vFast*1.60e6;%in Hz
            SweepFrequencyEnd = -vSlow*1.6e6;
            
            SweepReturnTime= SweepTimeEnd +1e-3;%return the frequency to the starting point.
            SweepRate = (SweepFrequencyEnd-SweepFrequencyStart)/(SweepTimeEnd - SweepTimeStart);
            ReturnSweepRate = (SweepFrequencyStart-SweepFrequencyEnd)/(SweepReturnTime - SweepTimeEnd);
            FluorvTime = zeros(4,round((SweepTimeEnd-SweepTimeStart/hFine))+round((SimTime-SweepTimeEnd-SweepTimeStart)/hCoarse));
            CubeFluorvTime = zeros(2,round((SweepTimeEnd-SweepTimeStart/hFine))+round((SimTime-SweepTimeEnd-SweepTimeStart)/hCoarse));
            while time <= SimTime
                vStart = XV(6,:);
                
                if (time-(SweepTimeStart-hCoarse)>=0 && (SweepTimeEnd+hCoarse-time)>0)
                    %We are slowing.
                    h = hFine;
                    %Symplectic Integrator
                    acc = LongBeamSlowingAcc(XV,time,SweepTimeStart,SweepRate,...
                        SweepFrequencyEnd,lambda0,lambda1,SweepTimeEnd,SweepFrequencyStart,ReturnSweepRate,SweepReturnTime);
                    XV(6,:) = XV(6,:) + acc.*h/2;
                    XV(1:3,:) = XV(1:3,:) + XV(4:6,:) *h;
                    acc = LongBeamSlowingAcc(XV,time+h/2,SweepTimeStart,SweepRate,...
                        SweepFrequencyEnd,lambda0,lambda1,SweepTimeEnd,SweepFrequencyStart,ReturnSweepRate,SweepReturnTime);
                    XV(6,:) = XV(6,:) + acc.*h/2;
                    
                else
                    %We are not slowing.
                    h = hCoarse;
                    %acc = ElectricGuideAcc(XV,mass);
                    NumberOfMolecules = length(XV(9,:));
                    acc = zeros(2,NumberOfMolecules);
                    XV(4:5,:) = XV(4:5,:) + acc.*h/2;
                    XV(1:3,:) = XV(1:3,:) + XV(4:6,:) *h;
                    %acc = ElectricGuideAcc(XV,mass);
                    %         acc = zeros(2,NumberOfMolecules);
                    %         XV(4:5,:) = XV(4:5,:) + acc.*h/2;
                    
                    XV(10,:) =XV(10,:) +((XV(1,:).^2+XV(2,:).^2<0.005^2)&(XV(1,:).^2+(XV(3,:)-0.4).^2<0.005^2)&...
                        (XV(2,:).^2+(XV(3,:)-0.4).^2<0.005^2)&XV(4,:)<=10&XV(5,:)<=10&XV(6,:)<=10);  %Have an index which increases for each molecule that is trappable.0 Not trapped.  1 trapped 2 was trapped (don't double count.
                    NumberOfTrappable =NumberOfTrappable+sum( XV(10,:)==1) ;
                end
                
                StateStart = XV(9,:);
                
                [XV, LostToDelta] = LossIntoDeltaAndv2(XV,vStart,StateStart,v2angle,v2BeamSize,LostToDelta);
                
                SimIndex = SimIndex + 1;
                time = time+h;
                
                %waitbar(time/SimTime)
%                 
%                 if mod(SimIndex,500) == 1
%                     figure(1)
%                     subplot(2,1,1)
%                     plot(XV(3,:),XV(2,:),'.')
%                     axis([0 0.42 -0.05 0.05])
%                     xlabel('z [m]')
%                     ylabel('x [m]')
%                     subplot(2,1,2)
%                     plot(XV(6,:),XV(8,:),'.',[0 150],[0 150],'-r')
%                     axis([0 120 -0.05 120])
%                     xlabel('Current v [m/s]')
%                     ylabel('Initial Velocity [m/s]')
%                     NumberOfTrappable
%                     time
%                 end
                
                XV = XV(:,XV(3,:)<CellToMotDist+.1 &abs(XV(1,:))<MOTsize*2.5 & abs(XV(2,:))<MOTsize*2.5); %Take out molecules that are too far away to make a difference.
                
                
                
                FluorvTime(1,SimIndex) = time; %time, motregion, trappable, number in v=2
%                 CubeFluorvTime(1,SimIndex) = time;
%                 DetectableCube = zeros(1,size(XV,2)) + ((XV(3,:)-CellToCubeDist-tan(detectionAngle)*XV(1,:)).^2+(XV(2,:).^2) <= CubeBeamSize^2 & XV(9,:) == 1); %In the detection region and in v = 0 or 1
%                 CountsCube = zeros(1,size(XV,2)) + round(h*(s0/(1+s0))*gammaV0./(2*(2+4/(gammaV0^2*(1+s0))*(motDetuning(l)+f0*XV(6,:)*sin(detectionAngle)/c).^2))).*DetectableCube;
%                 CubeFluorvTime(2,SimIndex) = sum(CountsCube);
%                 
                DetectableMOT = zeros(1,size(XV,2)) + ((XV(3,:)-CellToMotDist-tan(detectionAngle)*XV(1,:)).^2+(XV(2,:).^2) <= MOTsize^2 & XV(9,:) == 1); %In the detection region and in v = 0 or 1
%                 CountsMOT = zeros(1,size(XV,2)) + round(h*(s0/(1+s0))*gammaV0./(2*(2+4/(gammaV0^2*(1+s0))*(motDetuning(l)+f0*XV(6,:)*sin(detectionAngle)/c).^2))).*DetectableMOT;
%                CountsMOT = zeros(1,size(XV,2)) +round(h*(s0*gammaV0/2)./(1+s0+(4/gammaV0^2)*((motDetuning(l)+(f0+motDetuning(l))*sqrt((1+XV(6,:)*sin(detectionAngle)/c)/(1-XV(6,:)*sin(detectionAngle)/c)).^2)))).*DetectableMOT;
               CountsMOT = zeros(1,size(XV,2)) +round(h*(s0*gammaV0/2)./(1+s0 +(4/gammaV0^2)*(motDetuning(l)+(f0+motDetuning(l))*XV(6,:)*sin(detectionAngle)/c).^2)).*DetectableMOT;
               
               FluorvTime(2,SimIndex) = sum(CountsMOT);
                %                 FluorvTime(2,SimIndex) = length(find((XV(2,:).^2+XV(1,:).^2+(XV(3,:)-CellToMotDist).^2)<MOTsize^2 & XV(9,:)==1));  %Assuming v =2 repump is off in the MOT detection.
                %         FluorvTime(2,SimIndex) = length(find((XV(2,:).^2+XV(1,:).^2+(XV(3,:)-CellToMotDist).^2)<MOTsize^2)); %assume v=2 on in MOT region
                
                %    FluorvTime(2,SimIndex) = length(find(sqrt(XV(2,:).^2+XV(1,:).^2)<0.005...
                %       & abs(XV(3,:)-0.4)<0.005 ));  %Assuming v =2 repump is on in the MOT detection.
                
                FluorvTime(3,SimIndex) = NumberOfTrappable;
                v2 = find(XV(:,XV(9,:)==0));
                FluorvTime(4,SimIndex) = sum(v2);
                
                
            end
            LostToDelta
            NumberOfTrappable
            Trapped(i) =NumberOfTrappable;
   
%             figure()

%             plot(FluorvTime(1,:),FluorvTime(2,:))
%            
%             xlabel('time [s]')
%             xlim([0 .025]);
%             ylabel('Fluorescence [a.u.]')
%             plot(FluorvTime(1,:),FluorvTime(2,:))
%             plot(FluorvTime(1,:),FluorvTime(4,:))
%             xlabel('time [s]')
%             xlim([0 .025]);
%             ylabel('Number in v=2')
% 
%             title([num2str(vF(i)) ' m/s with spread of ' num2str(vSpread(j)) ', detection angle at ' num2str(detectionAngle) ', laser detuning ' num2str(motDetuning(l)/10^6) ' MHz (detect ' num2str(detectVels(l)) 'm/s'])
%            
%             save([num2str(k) '_vF' num2str(vF(i)) '_' num2str(vSpread(j)*10)], 'FluorvTime', 'CubeFluorvTime');
            save([num2str(k) '_vF' num2str(vF(i)) '_' num2str(vSpread(j)*10) '_' num2str(detectVels(l))], 'FluorvTime');
            %            
        end
    end
end
end
%filename = ['Sept25FvTStart' num2str(SweepTimeStart*10^3) 'msChirpLgth' ...
%     num2str((SweepTimeEnd-SweepTimeStart)*10^3) 'msFrom' num2str(vFast) ...
%     'to' num2str(vSlow) 'v2ClosedLasersNormalScan']
% %filename = 'Sept25v2Open'
%save('SlowingResults\Sept19NoAcc','FluorvTime')

%save(['SlowingResults\' filename],'FluorvTime')
