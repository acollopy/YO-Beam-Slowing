function [detuningTime] = GenerateFrequencySweep(vFast, vSlow, sweepTimeStart, sweepTimeEnd, hFine, hCoarse, simTime)
%GenerateFrequencySweep initializes the time vector and the frequency
%detuning (for v =1,2) at each point in time for the slowing sequence

sweepFrequencyStart = -vFast*1.60e6;%in Hz
sweepFrequencyEnd = -vSlow*1.6e6;

detuningTime(1,:) = [ 0:hCoarse:sweepTimeStart (sweepTimeStart + hFine):hFine:sweepTimeEnd (sweepTimeEnd + hCoarse):hCoarse:simTime];
detuningTime(2,:) = [ sweepFrequencyStart*ones(1,sweepTimeStart/hCoarse +1) sweepFrequencyStart+(sweepFrequencyEnd-sweepFrequencyStart)/((sweepTimeEnd-sweepTimeStart)/hFine):(sweepFrequencyEnd-sweepFrequencyStart)/((sweepTimeEnd-sweepTimeStart)/hFine):sweepFrequencyEnd sweepFrequencyEnd*ones(1,(simTime-sweepTimeEnd)/hCoarse)];
% figure
% plot(detuningTime(1,:), [1:size(detuningTime,2)])
% plot(detuningTime(1,:), detuningTime(2,:))

end

