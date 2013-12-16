function [acc] = LongBeamSlowingAcc(XV,time,SweepTimeStart,SweepRate,SweepFrequencyEnd,lambda0,lambda1,SweepTimeEnd,...
    SweepFrequencyStart,ReturnSweepRate,SweepReturnTime)

    BeamSize = .005;
    DeltaLaser1 =  (SweepFrequencyStart+(time-SweepTimeStart).*SweepRate).*sign(hvs(time-SweepTimeStart).*hvs(SweepTimeEnd-time))...
       +(SweepFrequencyEnd+(time-SweepTimeEnd).*ReturnSweepRate).*sign(hvs(time-SweepTimeEnd).*hvs(SweepReturnTime-time));
    DeltaDoppler1 = 2*pi.*XV(6,:).*(1/lambda1)+2*pi*DeltaLaser1;
%     AccInd = (DeltaDoppler+2*pi*200e6)./(2*pi*1000);

    DeltaLaser0 = DeltaLaser1*648/614;
    DeltaDoppler0 = 2*pi.*XV(6,:).*(1/lambda0)+2*pi*DeltaLaser0;
    %AccFullPower = -3.18e13.*(4.70e14+DeltaDoppler.^2).^(-1)*194604; %Calculated in matehamtica BeamSlowing3LevelSystem.nb
    AccFullPower = -1.09e14.*(2.206e15+DeltaDoppler0.^2+4.84*DeltaDoppler1.^2).^(-1)*194604;
    if (time-SweepTimeStart >0 && SweepTimeEnd>time)
        acc = AccFullPower.*XV(9,:).*hvs(XV(3,:)).*hvs(BeamSize-sqrt(XV(1,:).^2+XV(2,:).^2));
    else
        acc = 0;
end