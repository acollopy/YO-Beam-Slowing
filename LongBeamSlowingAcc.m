function [acc] = LongBeamSlowingAcc(XV,time,SweepTimeStart,SweepFrequencyEnd,lambda0,lambda1,SweepTimeEnd,...
    SweepFrequencyStart,returnSweepLength,longBeamSize)

SweepReturnTime= SweepTimeEnd +returnSweepLength;%return the frequency to the starting point.
SweepRate = (SweepFrequencyEnd-SweepFrequencyStart)/(SweepTimeEnd - SweepTimeStart);
ReturnSweepRate = (SweepFrequencyStart-SweepFrequencyEnd)/(returnSweepLength);

DeltaLaser1 =  (SweepFrequencyStart+(time-SweepTimeStart).*SweepRate).*ceil(time-SweepTimeStart).*ceil(SweepTimeEnd-time)...
    +(SweepFrequencyEnd+(time-SweepTimeEnd).*ReturnSweepRate).*ceil(time-SweepTimeEnd).*ceil(SweepReturnTime-time);
DeltaDoppler1 = 2*pi.*XV(6,:).*(1/lambda1)+2*pi*DeltaLaser1;
%     AccInd = (DeltaDoppler+2*pi*200e6)./(2*pi*1000);

DeltaLaser0 = DeltaLaser1*648/614;
DeltaDoppler0 = 2*pi.*XV(6,:).*(1/lambda0)+2*pi*DeltaLaser0;
%AccFullPower = -3.18e13.*(4.70e14+DeltaDoppler.^2).^(-1)*194604; %Calculated in matehamtica BeamSlowing3LevelSystem.nb
AccFullPower = -1.09e14.*(2.206e15+DeltaDoppler0.^2+4.84*DeltaDoppler1.^2).^(-1)*194604;
if (time-SweepTimeStart >0 && SweepTimeEnd>time)
    NumberOfMolecules = size(XV,2);
    InSlowingBeam = zeros(1,NumberOfMolecules)+(XV(3,:)>0 & XV(1,:).^2+XV(2,:).^2 <= longBeamSize^2 & XV(9,:) == 1); %in beam and in state 1
    acc = AccFullPower.*XV(9,:).*InSlowingBeam;
    %         acc = AccFullPower.*XV(9,:).*hvs(XV(3,:)).*hvs(BeamSize-sqrt(XV(1,:).^2+XV(2,:).^2));
else
    acc = 0;
end
end