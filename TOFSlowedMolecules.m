clear

%setup constants
gammaV0 = 2*pi*5*10^6;
kB = 1.38e-23;
mass = 105*1.6605e-27;
T = 3.2;
CellToMotDist = 0.44;
CellToCubeDist = .22;
CubeBeamSize = .002;
MOTsize = .010;
longBeamSize = .01;
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
returnSweepLength = .001; %ramp back to first frequency in 1 ms


% motDetuning = 0*10^6; %Hz, compared to absolute frequency
% detectVels = [50 60 70 80 90 100 110 120 130];
detectVels = [40 60];
%motDetuning = -1.6*detectVels*10^6;

motDetuning = -f0*detectVels*sin(detectionAngle)/c;
averagenum = [1 1];

vF = [80];
vSpread = [16];
% vSpread = sqrt(kB*T/mass);

v2BeamSize = .005; %radius

hFine = 1e-6; %simulation time step
hCoarse = 1e-4;
simTime = 40e-3;

                vFast = 150;
                vSlow = 150;
                sweepTimeStart =25e-3;
                sweepTimeEnd = 25.01e-3;
                %             vSlow = vFast(j) - deltaV(k);
                
                SweepFrequencyStart = -vFast*1.60e6;%in Hz
                SweepFrequencyEnd = -vSlow*1.6e6;

detuningTime = GenerateFrequencySweep(vFast,vSlow,sweepTimeStart,sweepTimeEnd,hFine,hCoarse,simTime); %first row is time, second is detuning from resonance for v=1,2 (Hz)

for l = 1:size(detectVels,2) %go through velocity detunings
    for k = averagenum(1):averagenum(2) % go through process multiple times for averaging purposes
        for i =1:size(vF,2) % go through different forward velocities
            for j = 1:size(vSpread,2) % go through different forward velocity spreads
                rng('shuffle');

                XV = SetUpInitDistro(vF(i),vSpread(j));
                NumberOfMolecules = size(XV,2);

                %Frequency Sweep Parameters.
                


                vint = zeros(3,NumberOfMolecules);
                acc = zeros(3,NumberOfMolecules);
                SimIndex = 1;
                FluorvTime = zeros(4,round((sweepTimeEnd-SweepTimeStart/hFine))+round((simTime-sweepTimeEnd-SweepTimeStart)/hCoarse));
%                 CubeFluorvTime = zeros(2,round((SweepTimeEnd-SweepTimeStart/hFine))+round((SimTime-SweepTimeEnd-SweepTimeStart)/hCoarse));
               

                NumberOfMolecules=size(XV,2);
                vStart = zeros(1,NumberOfMolecules);
                NumberOfTrappable = 0;
                LostToDelta = 0;
                

                
                for t = 1:size(detuningTime)
                    vStart = XV(6,:);
                    
                    if (time-(SweepTimeStart-hCoarse)>=0 && (sweepTimeEnd+hCoarse-time)>0)
                        %We are slowing.
                        h = hFine;
                        %Symplectic Integrator
                        acc = LongBeamSlowingAcc(XV,time,SweepTimeStart,...
                            SweepFrequencyEnd,lambda0,lambda1,sweepTimeEnd,SweepFrequencyStart,returnSweepLength,longBeamSize);
                        XV(6,:) = XV(6,:) + acc.*h/2;
                        XV(1:3,:) = XV(1:3,:) + XV(4:6,:) *h;
                        acc = LongBeamSlowingAcc(XV,time+h/2,SweepTimeStart,...
                            SweepFrequencyEnd,lambda0,lambda1,sweepTimeEnd,SweepFrequencyStart,returnSweepLength,longBeamSize);
                        XV(6,:) = XV(6,:) + acc.*h/2;
                        
                    else
                        %We are not slowing.
                        h = hCoarse;
                        %acc = ElectricGuideAcc(XV,mass);
                        NumberOfMolecules = size(XV,2);
                        acc = zeros(2,NumberOfMolecules);
                        XV(4:5,:) = XV(4:5,:) + acc.*h/2;
                        XV(1:3,:) = XV(1:3,:) + XV(4:6,:) *h;
                        
                        XV(10,:) =XV(10,:) +((XV(1,:).^2+XV(2,:).^2<0.005^2)&(XV(1,:).^2+(XV(3,:)-0.4).^2<0.005^2)&...
                            (XV(2,:).^2+(XV(3,:)-0.4).^2<0.005^2)&XV(4,:)<=10&XV(5,:)<=10&XV(6,:)<=10);  %Have an index which increases for each molecule that is trappable.0 Not trapped.  1 trapped 2 was trapped (don't double count.
                        NumberOfTrappable =NumberOfTrappable+sum( XV(10,:)==1) ;
                    end
                   
                    
                    [XV, LostToDelta] = LossIntoDeltaAndv2(XV,vStart,v2angle,v2BeamSize,LostToDelta);
                    
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

                    DetectableMOT = zeros(1,size(XV,2)) + ((XV(3,:)-CellToMotDist-tan(detectionAngle)*XV(1,:)).^2+(XV(2,:).^2) <= MOTsize^2 & XV(9,:) == 1); %In the detection region and in v = 0 or 1
                    CountsMOT = zeros(1,size(XV,2)) +round(h*(s0*gammaV0/2)./(1+s0 +(4/gammaV0^2)*(motDetuning(l)+(f0+motDetuning(l))*XV(6,:)*sin(detectionAngle)/c).^2)).*DetectableMOT;
                    
                    FluorvTime(2,SimIndex) = sum(CountsMOT);
                    
                    FluorvTime(3,SimIndex) = NumberOfTrappable;
                    v2 = find(XV(:,XV(9,:)==0));
                    
                    FluorvTime(4,SimIndex) = sum(v2);
                    
                    
                end
                LostToDelta
                NumberOfTrappable
                Trapped(i) =NumberOfTrappable;
                
                            figure()
                
                            plot(FluorvTime(1,:),FluorvTime(2,:))
                
                            xlabel('time [s]')
                            xlim([0 .025]);
                            ylabel('Fluorescence [a.u.]')

                            title([num2str(vF(i)) ' m/s with spread of ' num2str(vSpread(j)) ', detection angle at ' num2str(detectionAngle) ', laser detuning ' num2str(motDetuning(l)/10^6) ' MHz (detect ' num2str(detectVels(l)) 'm/s'])
               save([num2str(k) '_vF' num2str(vF(i)) '_' num2str(vSpread(j)*10) '_' num2str(detectVels(l))], 'FluorvTime');
                %
            end
        end
    end
end
