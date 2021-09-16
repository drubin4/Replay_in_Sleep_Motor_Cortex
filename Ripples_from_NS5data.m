%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021


%%This script is going to look overnight at Ripples, and look
%%for a relationship to replay events. 

%This will also incorporate the NS5 data that was separately compiled

%% Open intracranial neural data as well as surface EEG data:

clear all
close all
clc

%% Load the kinematic data

tic
load('AnalysisWorkspace.mat')
load('ProcessedKalmanModelData.mat')
load('STCEs_18Speeds_21Thresholds.mat')
toc
%% Load the NS5 data

rippledir = 'Ripples';
tic
load([rippledir '/' 'LateralArrayAllRipples_lite.mat']);
toc
tic
load([rippledir '/' 'MedialArrayAllRipples_lite.mat']);
toc

[MedialTimeAxis, ind] = sort(MedialTimeAxis);
MedialArrayRippleYesNo = MedialArrayRippleYesNo(:,ind);

[LateralTimeAxis, ind] = sort(LateralTimeAxis);
LateralArrayRippleYesNo = LateralArrayRippleYesNo(:,ind);

%%
makefigs = 0;       %Do (1) or do not (0) print figures as the analysis progresses

%% Make the sleep metric:

clear VeryInterestingIndices NumberVeryInterestIndices durationofcomparisonperiod totalfiringrateduringperiod

EEGSpectrogramTimeAxis = datenum(EEG_Spectrogram_DateTimeAxis);
ThresholdFreq = zeros(1,length(EEGSpectrogramTimeAxis));
ThresholdPercent = 95;
for t = 1:length(EEGSpectrogramTimeAxis)
    
    ThresholdProp = ThresholdPercent/100;
    LowPassThreshold = find(EEG_Spectrogram_Frequency_Axis ==20);
    spectrogram_slice = EEG_Spectrogram(1:LowPassThreshold,t);
    spectrogram_slice = exp(spectrogram_slice);
    spectrogram_slice_norm = spectrogram_slice./(sum(spectrogram_slice));
    integer_spectrogram_slice_norm = cumsum(spectrogram_slice_norm);
    cutoffInd = find(integer_spectrogram_slice_norm<ThresholdProp);
    cutoffInd = cutoffInd(end);
    
    ThresholdFreq(t) = EEG_Spectrogram_Frequency_Axis(cutoffInd);
    
end
smoothedthresholdFrequency = movmean(ThresholdFreq,[60 0]);  %create 30-second moving average

%% Find times of each replay event at each speed:
    
speeds = [6 12 13 14]; 
CCrossingsX = zeros(4,length(DateNum_IC_TimeAxis));

for z = 1:4
    TheseAreTheIndicesOfthe_STCEs = CCAllInterestingIndicesModified{1,speeds(z)};
    
    CCrossings = DateNum_IC_TimeAxis*0;
    CCrossings(TheseAreTheIndicesOfthe_STCEs) = 2;
    CCrossings = CCrossings+z*2;
    
    CCrossingsX(z,:) = CCrossings;
end

%% Plot some examples with EEG, spike Raster, kalman filter and NS5 ripple activity during some replay events

i = 11; j = 6; k = 77;

IndicesToLookAt = CCAllInterestingIndicesModified{i,j};

STI_IC = IndicesToLookAt(k) - 50*1; 
ETI_IC = IndicesToLookAt(k) + 50*10;

% Make the time segment
  
Xlim1 = TimeAxis(STI_IC);
Xlim2 = TimeAxis(ETI_IC);

Xlim1t = TimeAxis(STI_IC)-minutes(0.25);
Xlim2t = TimeAxis(ETI_IC)+minutes(0.25);

StartTimeIndex_ICT = find(TimeAxis > Xlim1t);
StartTimeIndex_ICT = StartTimeIndex_ICT(1);

EndTimeIndex_ICT = find(TimeAxis > Xlim2t);
EndTimeIndex_ICT = EndTimeIndex_ICT(1)-1;

% Plot the EEG data for a minute around the ripple
Xlim1_eeg = TimeAxis(StartTimeIndex_ICT-(50*60));
Xlim2_eeg = TimeAxis(EndTimeIndex_ICT+(50*60));

% Find the Xlim parameters for the Ripple NS5 data

StartTimeIndex_MedialArrayNS5 = find(MedialTimeAxis > Xlim1t);
StartTimeIndex_MedialArrayNS5 = StartTimeIndex_MedialArrayNS5(1);

EndTimeIndex_MedialArrayNS5 = find(MedialTimeAxis > Xlim2t);
EndTimeIndex_MedialArrayNS5 = EndTimeIndex_MedialArrayNS5(1)-1;

StartTimeIndex_LateralArrayNS5 = find(LateralTimeAxis > Xlim1t);
StartTimeIndex_LateralArrayNS5 = StartTimeIndex_LateralArrayNS5(1);

EndTimeIndex_LateralArrayNS5 = find(LateralTimeAxis > Xlim2t);
EndTimeIndex_LateralArrayNS5 = EndTimeIndex_LateralArrayNS5(1)-1;
%

fignum = 131;

h = figure(fignum);
clf
set(fignum,'Position',[1 -145 1680 933])
set(fignum,'DefaultAxesFontSize',16);
h.Color = [1 1 1];

subplot(5,1,1)
h = plotEEGSpectrogram(Xlim1_eeg,Xlim2_eeg,fignum,EEGSpectrogramTimeAxis,EEG_Spectrogram_Frequency_Axis,EEG_Spectrogram);
title('Scalp EEG Spectrogram')

subplot(5,1,2)
S = plotRasterplot(StartTimeIndex_ICT,EndTimeIndex_ICT,fignum,DateNum_IC_TimeAxis,Thresholded_zSP);
title('Spike Raster')

h = subplot(5,1,3);
plot(DateNum_IC_TimeAxis(StartTimeIndex_ICT:EndTimeIndex_ICT),Kalman_Data(1:2,(StartTimeIndex_ICT):(EndTimeIndex_ICT)))
axis tight; box off
datetick('x','keeplimits');
title('Kalman Filter Output')



subplot(5,1,4)
hold on
S = plotRippleRaster(StartTimeIndex_MedialArrayNS5,EndTimeIndex_MedialArrayNS5,fignum,MedialTimeAxis,MedialArrayRippleYesNo');
title('Medial Array Ripple')

subplot(5,1,5)
hold on
S = plotRippleRaster(StartTimeIndex_LateralArrayNS5,EndTimeIndex_LateralArrayNS5,fignum,LateralTimeAxis,LateralArrayRippleYesNo');
title('Lateral Array Ripple')


PrintFigTitle = ['Example With Ripple Presence.svg'];
%print(PrintFigTitle,'-dsvg','-painters')


%% Overnight, how much of night is in Ripple?
% Define "Ripple" by the sum of the ripple count on each array
clear percentcorr1x percentcorr2x percentcorr3x percentcorr4x PercentRippleTime
RipCounter = 1;

%%% Downsample all the NS5 data again:

MedialArrayRippleYesNo_ds = downsample(MedialArrayRippleYesNo',10);
LateralArrayRippleYesNo_ds = downsample(LateralArrayRippleYesNo',10);
MedialArrayRippleYesNo_ds = MedialArrayRippleYesNo_ds';
LateralArrayRippleYesNo_ds = LateralArrayRippleYesNo_ds';
MedialTimeAxis_ds = downsample(MedialTimeAxis,10);
LateralTimeAxis_ds = downsample(LateralTimeAxis,10);

%Need to create an indexing for TimeAxis (50Hz) into MedialTimeAxis (1000Hz) and
%LateralTimeAxis(1000Hz)
%numericTimeAxis = datenum(TimeAxis);
%numericMedialTimeAxis = datenum(MedialTimeAxis_ds);
%numericLateralTimeAxis = datenum(LateralTimeAxis_ds);

TimeAxis_to_MedialTimeAxis = zeros(size(TimeAxis));
TimeAxis_to_LateralTimeAxis = zeros(size(TimeAxis));
%%
'starting'

tic

MedialTimeAxis_ds_adj = MedialTimeAxis_ds - seconds(0.7);

for i = 1  %length(numericTimeAxis)
    
    [d mind] = min(abs(MedialTimeAxis_ds_adj - TimeAxis(i)));
    TimeAxis_to_MedialTimeAxis(i) = mind;
    
  %  [d mind] = min(abs(numericLateralTimeAxis - numericTimeAxis(i)));
  %  numericTimeAxis_to_LateralTimeAxis(i) = mind;
    
end

for i = 2:length(TimeAxis)
    
    if abs(seconds(TimeAxis(i-1) - MedialTimeAxis_ds_adj(TimeAxis_to_MedialTimeAxis(i-1))) < 1)
        
        try
        [d mind] = min(abs(MedialTimeAxis_ds_adj(TimeAxis_to_MedialTimeAxis(i-1):(TimeAxis_to_MedialTimeAxis(i-1)+200)) - TimeAxis(i)));
        TimeAxis_to_MedialTimeAxis(i) = TimeAxis_to_MedialTimeAxis(i-1)+mind-1;
        catch 
        [d mind] = min(abs(MedialTimeAxis_ds_adj - TimeAxis(i)));
        TimeAxis_to_MedialTimeAxis(i) = mind;
        end    
        
    
    else
        [d mind] = min(abs(MedialTimeAxis_ds_adj - TimeAxis(i)));
        TimeAxis_to_MedialTimeAxis(i) = mind;
    
    end

    
    
  %  [d mind] = min(abs(numericLateralTimeAxis - numericTimeAxis(i)));
  %  numericTimeAxis_to_LateralTimeAxis(i) = mind;
    
end
toc

%save('TimeAxes_downsampled_adj.mat','TimeAxis_to_MedialTimeAxis','-v7.3')
toc
%%

clear PercentRippleTime percentcorr1x percentcorr2x percentcorr3x percentcorr4x
RipCounter = 1;

for RippleThreshold = [1 2 4 8 16]

Xlim1t =  datetime(2020,11,24,23,0,0);
Xlim2t =  datetime(2020,11,25,9,0,0);

StartTimeIndex_ICT = find(TimeAxis > Xlim1t);
StartTimeIndex_ICT = StartTimeIndex_ICT(1);
EndTimeIndex_ICT = find(TimeAxis > Xlim2t);
EndTimeIndex_ICT = EndTimeIndex_ICT(1)-1;
Xlim1_eeg = TimeAxis(StartTimeIndex_ICT);
Xlim2_eeg = TimeAxis(EndTimeIndex_ICT);

% Find the Xlim parameters for the Ripple NS5 data

StartTimeIndex_MedialArrayNS5 = find(MedialTimeAxis_ds_adj > Xlim1t);
StartTimeIndex_MedialArrayNS5 = StartTimeIndex_MedialArrayNS5(1);

EndTimeIndex_MedialArrayNS5 = find(MedialTimeAxis_ds_adj > Xlim2t);
EndTimeIndex_MedialArrayNS5 = EndTimeIndex_MedialArrayNS5(1)-1;
% 
% StartTimeIndex_LateralArrayNS5 = find(LateralTimeAxis > Xlim1t);
% StartTimeIndex_LateralArrayNS5 = StartTimeIndex_LateralArrayNS5(1);
% 
% EndTimeIndex_LateralArrayNS5 = find(LateralTimeAxis > Xlim2t);
% EndTimeIndex_LateralArrayNS5 = EndTimeIndex_LateralArrayNS5(1)-1;

%
AggregateMedialArrayRipple = (sum(MedialArrayRippleYesNo_ds));
%AggregateLateralArrayRipple = (sum(LateralArrayRippleYesNo_ds));

RippleHere_M = AggregateMedialArrayRipple>RippleThreshold;
RippleHere_M = double(RippleHere_M);


RippleTime = RippleHere_M(StartTimeIndex_MedialArrayNS5:EndTimeIndex_MedialArrayNS5);

PercentRippleTime(RipCounter) = sum(RippleTime)./numel(RippleTime);
clear Tdelay
counter=1;

    for z = 1:1:701
        tic
    OneXCrossings = (CCrossingsX(1,StartTimeIndex_ICT:EndTimeIndex_ICT) - 2)/2;
    TwoXCrossings = (CCrossingsX(2,StartTimeIndex_ICT:EndTimeIndex_ICT) - 4)/2;
    ThrXCrossings = (CCrossingsX(3,StartTimeIndex_ICT:EndTimeIndex_ICT) - 6)/2;
    FouXCrossings = (CCrossingsX(4,StartTimeIndex_ICT:EndTimeIndex_ICT) - 8)/2;

    OneXCrossings = circshift(OneXCrossings,z-351);
    TwoXCrossings = circshift(TwoXCrossings,z-351);
    ThrXCrossings = circshift(ThrXCrossings,z-351);
    FouXCrossings = circshift(FouXCrossings,z-351);

    Tdelay(counter) = z-351;

    a = (find(OneXCrossings == 1)) + StartTimeIndex_ICT;   
    percentcorr1x(RipCounter,counter) = sum(RippleHere_M(TimeAxis_to_MedialTimeAxis(a)))/length(a);

    a =  (find(TwoXCrossings == 1)) + StartTimeIndex_ICT;   
    percentcorr2x(RipCounter,counter) = sum(RippleHere_M(TimeAxis_to_MedialTimeAxis(a)))/length(a);

    a =  (find(ThrXCrossings == 1)) + StartTimeIndex_ICT; 
    percentcorr3x(RipCounter,counter) = sum(RippleHere_M(TimeAxis_to_MedialTimeAxis(a)))/length(a);
    
    a =  (find(FouXCrossings == 1)) + StartTimeIndex_ICT; 
    percentcorr4x(RipCounter,counter) = sum(RippleHere_M(TimeAxis_to_MedialTimeAxis(a)))/length(a);


    toc
    counter = counter+1
    end

RipCounter = RipCounter+1;

end


%%
TdelayS = Tdelay*0.02;

fignum = 2;
h=figure(fignum);
set(fignum,'Position',[1 33 1536 755])
set(fignum,'DefaultAxesFontSize',16);
h.Color = [1 1 1];

subplot(2,2,1)
hold on
plot(TdelayS,percentcorr1x)
a = get(gca,'colororder');
plot(TdelayS,PercentRippleTime(1)+0*TdelayS,'linestyle','--','Color',a(1,:))
plot(TdelayS,PercentRippleTime(2)+0*TdelayS,'linestyle','--','Color',a(2,:))
plot(TdelayS,PercentRippleTime(3)+0*TdelayS,'linestyle','--','Color',a(3,:))
plot(TdelayS,PercentRippleTime(4)+0*TdelayS,'linestyle','--','Color',a(4,:))
plot(TdelayS,PercentRippleTime(5)+0*TdelayS,'linestyle','--','Color',a(5,:))
%line([2 2],[0 1],'LineWidth',1,'Color',[0 0 0],'linestyle','--')
%plot(TdelayS,PercentRippleTime(6)+0*TdelayS,'linestyle','--','Color',a(6,:))
%plot(TdelayS,PercentRippleTime(7)+0*TdelayS,'linestyle','--','Color',a(7,:))
box off
axis tight
ylim([0 1])
xlim([-7 7])
legend('Ripple: >1 channel','Ripple: >2 channel','Ripple: >4 channel','Ripple: >8 channel','Ripple: >16 channel','Location','NorthWest')
xlabel('STCE Template Offset (seconds)')
ylabel('Proportion')
title('1x Replay Events Occuring During Ripple')

subplot(2,2,2)
hold on
plot(TdelayS,percentcorr2x)
a = get(gca,'colororder');
plot(TdelayS,PercentRippleTime(1)+0*TdelayS,'linestyle','--','Color',a(1,:))
plot(TdelayS,PercentRippleTime(2)+0*TdelayS,'linestyle','--','Color',a(2,:))
plot(TdelayS,PercentRippleTime(3)+0*TdelayS,'linestyle','--','Color',a(3,:))
plot(TdelayS,PercentRippleTime(4)+0*TdelayS,'linestyle','--','Color',a(4,:))
plot(TdelayS,PercentRippleTime(5)+0*TdelayS,'linestyle','--','Color',a(5,:))
%line([1 1],[0 1],'LineWidth',1,'Color',[0 0 0],'linestyle','--')
%plot(TdelayS,PercentRippleTime(6)+0*TdelayS,'linestyle','--','Color',a(6,:))
%plot(TdelayS,PercentRippleTime(7)+0*TdelayS,'linestyle','--','Color',a(7,:))
box off
axis tight
ylim([0 1])
xlim([-7 7])
xlabel('STCE Template Offset (seconds)')
ylabel('Proportion')
title('2x Replay Events Occuring During Ripple')

subplot(2,2,3)
hold on
plot(TdelayS,percentcorr3x)
a = get(gca,'colororder');
plot(TdelayS,PercentRippleTime(1)+0*TdelayS,'linestyle','--','Color',a(1,:))
plot(TdelayS,PercentRippleTime(2)+0*TdelayS,'linestyle','--','Color',a(2,:))
plot(TdelayS,PercentRippleTime(3)+0*TdelayS,'linestyle','--','Color',a(3,:))
plot(TdelayS,PercentRippleTime(4)+0*TdelayS,'linestyle','--','Color',a(4,:))
plot(TdelayS,PercentRippleTime(5)+0*TdelayS,'linestyle','--','Color',a(5,:))
%line([2/3 2/3],[0 1],'LineWidth',1,'Color',[0 0 0],'linestyle','--')
%plot(TdelayS,PercentRippleTime(6)+0*TdelayS,'linestyle','--','Color',a(6,:))
%plot(TdelayS,PercentRippleTime(7)+0*TdelayS,'linestyle','--','Color',a(7,:))
box off
axis tight
ylim([0 1])
xlim([-7 7])
xlabel('STCE Template Offset (seconds)')
ylabel('Proportion')
title('3x Replay Events Occuring During Ripple')

subplot(2,2,4)
hold on
plot(TdelayS,percentcorr4x)
a = get(gca,'colororder');
plot(TdelayS,PercentRippleTime(1)+0*TdelayS,'linestyle','--','Color',a(1,:))
plot(TdelayS,PercentRippleTime(2)+0*TdelayS,'linestyle','--','Color',a(2,:))
plot(TdelayS,PercentRippleTime(3)+0*TdelayS,'linestyle','--','Color',a(3,:))
plot(TdelayS,PercentRippleTime(4)+0*TdelayS,'linestyle','--','Color',a(4,:))
plot(TdelayS,PercentRippleTime(5)+0*TdelayS,'linestyle','--','Color',a(5,:))
%line([2/4 2/4],[0 1],'LineWidth',1,'Color',[0 0 0],'linestyle','--')
%plot(TdelayS,PercentRippleTime(6)+0*TdelayS,'linestyle','--','Color',a(6,:))
%plot(TdelayS,PercentRippleTime(7)+0*TdelayS,'linestyle','--','Color',a(7,:))
box off
axis tight
ylim([0 1])
xlim([-7 7])
xlabel('STCE Template Offset (seconds)')
ylabel('Proportion')
title('4x Replay Events Occuring During Ripple')




PrintFigTitle = 'Relationship between NS5 ripples and replay events.svg';
print(PrintFigTitle,'-dsvg','-painters')

