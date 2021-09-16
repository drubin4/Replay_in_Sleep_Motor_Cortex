%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021

function [StartTimeIndex_IC EndTimeIndex_IC DateNum_IC_TimeAxis] = xLimtoIndex(xLim1, xLim2, TimeAxis)

StartTimeIndex_IC = find(TimeAxis > xLim1);
StartTimeIndex_IC = StartTimeIndex_IC(1);
EndTimeIndex_IC = find(TimeAxis > xLim2);
EndTimeIndex_IC = EndTimeIndex_IC(1);
DateNum_IC_TimeAxis = datenum(TimeAxis);

end

