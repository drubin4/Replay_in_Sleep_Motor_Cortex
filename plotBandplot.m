%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021


function S = plotBandplot(StartTimeIndex_IC,EndTimeIndex_IC,fignum,DateNum_IC_TimeAxis,BandPower,c_lim)
%PLOT Rasterplot of intracranial neural data; presume 96 channels

figure(fignum)
% 
% [a,b] = size(Thresholded_zSP(StartTimeIndex_IC:EndTimeIndex_IC,:)');
% CO = zeros(a,b,3);
% for i = 1:3
%    CO(:,:,i) = 1 - (Thresholded_zSP(StartTimeIndex_IC:EndTimeIndex_IC,:)'); 
% end

logBP = log2(BandPower(StartTimeIndex_IC:EndTimeIndex_IC,:)');
%logBP = (BandPower(StartTimeIndex_IC:EndTimeIndex_IC,:)');

S = imagesc(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),1:96,logBP); 

%caxis(c_lim);
%colorbar

% view(2)
% S.EdgeColor = 'none';
% S.LineStyle = 'none';
axis tight;
datetick('x','keeplimits');
ylabel('Channel');
%xlabel('Time');
%title('Spike Raster')


end

