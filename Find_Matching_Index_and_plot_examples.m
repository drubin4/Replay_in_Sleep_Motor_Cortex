%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021

%% This script calculates Matching Index between firing at different times and plots examples, using STCEs as the times to compare

%% Open intracranial neural data as well as surface EEG data:

clear all
close all
clc

%% Load all the data"

load('AnalysisWorkspace.mat')

load('ProcessedKalmanModelData.mat')

load('STCEs_18Speeds_21Thresholds.mat')

SigmaBand = load('Band2_Array_1.mat');
SigmaBand_a1 = SigmaBand.Band2_array1;

SigmaBand = load('Band2_Array_2.mat');
SigmaBand_a2 = SigmaBand.Band2_array2;


GammaBand = load('Band4_Array_1.mat');
GammaBand_a1 = GammaBand.Band4_array1;

GammaBand = load('Band4_Array_2.mat');
GammaBand_a2 = GammaBand.Band4_array2;

% 
Band1Range = [0 7.8125]
Band2Range = [11.7188 15.6250]
Band3Range = [19.5313 31.2500]
Band4Range = [39.0625 125.0000]
Band5Range = [128.9063 250.0000]

makefigs = 0;       %Do (1) or do not (0) print figures as the analysis progresses

%% Identify when all the replay events are occurring at each speed/threshold:

b = [0.1 0.25 0.5 0.75 0.9];
a = [1.0 1.1 1.2 1.4 1.6 1.8 2 3 4 5 6 8 10];
tCF = [b a];

i = 1        %i goes from 1 to 21; where 1 is the regular threshold and higher numbers are stricter threshold (higher correlation)
j = 6        %j goes from 1 to 18; j = 6 is 1x, j = 12 is 2x, j = 13 is 3x, j = 14 is 4x;


%TCEs is the index of each threshold crossing event at threshold level i,
%for speed tCF(j)

TCEs = CCAllInterestingIndicesModified{i,j};

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

%% Next find the average power in the sigma and gamma bands
% Sigma power will be for "spindles" and gamma power will be for "ripples"

meanSigma_a1 = mean(SigmaBand_a1,2);
meanSigma_a2 = mean(SigmaBand_a2,2);

meanGamma_a1 = mean(GammaBand_a1,2);
meanGamma_a2 = mean(GammaBand_a2,2);

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



%% Load the Spike Power and raw ncTX from medial array for analysis below looking at spike ordering:

DataSetToOpen1 = 'ncTX_Array_2.mat';
DataSetToOpen2 = 'SpikePower_Array_2.mat';

load(DataSetToOpen1)
load(DataSetToOpen2)

SP_nctx = single(ncTX_array2);                  %Array 2 is medial array
clear ncTX_array2
SP_sp = SpikePower_array2;
clear SpikePower_array2

% Z-score and threshold the IC data for plotting Rasters (these are also for display; no analysis is done on these data):
clear zSP_ncTX zSP_sp
%z-score each channel
for i = 1:96
    
   SPi = SP_nctx(:,i);  
   zSP_ncTX(:,i) = (SPi - mean(SPi))/std(SPi);   
    
   SPi = SP_sp(:,i);  
   zSP_sp(:,i) = (SPi - mean(SPi))/std(SPi);
   
end

%% Identify an order of unit acitivity based on spiking peak within each task-associated event, and then order the spiking activity during the replay events to see if it matches.
i = 1;
j = 6; %(for 1x speed)

IndicesToLookAt = CCAllInterestingIndicesModified{i,j};

%Game playing:
xLim1 = datetime(2020,11,24,14,50,0);
% xLim2 = datetime(2020,11,24,17,15,0);
% %Overnight    
% xLim1 = datetime(2020,11,24,23,0,0);
xLim2 = datetime(2020,11,25,9,0,0);
    
Indices_of_actualGamePlayEvents = TimeAxis(IndicesToLookAt);

a = find(Indices_of_actualGamePlayEvents>xLim1 & Indices_of_actualGamePlayEvents<xLim2);

Indices_of_actualGamePlayEvents = IndicesToLookAt(a);

ActualGamePlaySpikeOrder = zeros(length(Indices_of_actualGamePlayEvents),96);


sigmaF = .15;   
oSP = zSP_ncTX;
%oSP = zSP_sp;
    tic
    for j = 1:96
        oSP_j = oSP(:,j);
        smoothed_oSP_j = conv(oSP_j,exp(-(((-10:0.1:10).^2)./(2*sigmaF.^2))),'same');
        oSP(:,j) = smoothed_oSP_j;
    end
    toc
    
    
for k = 1:length(Indices_of_actualGamePlayEvents)
%for k = 50:70
STI_IC = Indices_of_actualGamePlayEvents(k) - 50*1; 
ETI_IC = Indices_of_actualGamePlayEvents(k) + 50*5;

%Find the order of firing based on time of max firing rate (try spikePower
%and ncTX separately and see which works better) within the X seconds of
%STCE
    
    TimeWindowToLookForPeak = 3;
    
    oSPi = oSP((Indices_of_actualGamePlayEvents(k)+25):(Indices_of_actualGamePlayEvents(k)+50*TimeWindowToLookForPeak),:);
    
    clear maxSPindex
    for j = 1:96
       oSP_j = oSPi(:,j);
       [a b] = max(oSP_j);
       b=b(1);
       maxSPindex(j) = b;
    end

    [a MaxFiringRateOrder] = sort(maxSPindex,'descend');
    
    [a FiringRateOrder] = sort(maxSPindex,'ascend');
    ActualGamePlaySpikeOrder(k,:) = FiringRateOrder;
    
%now plot the spiking data in order:
fignum = k*0+1;
figure(fignum)
set(fignum,'Position',[488 -134.2000 560 896.2000])
plotplot = 0;

if plotplot == 1
subplot(4,1,1)
OrderedThresholded_zSP = Thresholded_zSP(:,MaxFiringRateOrder);
S = plotRasterplot(STI_IC,ETI_IC,fignum,DateNum_IC_TimeAxis,OrderedThresholded_zSP);
%plot(DateNum_IC_TimeAxis(STI_IC:ETI_IC),mean(SP(STI_IC:ETI_IC,:),2))
t_lim = [datenum(DateNum_IC_TimeAxis(STI_IC)) datenum(DateNum_IC_TimeAxis(ETI_IC))];
xlim(t_lim);
datetick('x','keeplimits');
%ylabel('Mean Firing Rate')
pause(0.1)

subplot(4,1,2)
hold on
plot(DateNum_IC_TimeAxis(STI_IC:ETI_IC),Kalman_Data(1:2,STI_IC:ETI_IC))
plot([DateNum_IC_TimeAxis(Indices_of_actualGamePlayEvents(k)) DateNum_IC_TimeAxis(Indices_of_actualGamePlayEvents(k))],[-3 3],'g--')
axis tight; box off
datetick('x','keeplimits');
title('Kalman Filter Output')

subplot(4,1,3)
OrderedzSP_ncTX = zSP_ncTX(:,MaxFiringRateOrder);
S = plotRasterplot(STI_IC,ETI_IC,fignum,DateNum_IC_TimeAxis,OrderedzSP_ncTX);
%S = plotBandplot(STI_IC,ETI_IC,k,DateNum_IC_TimeAxis,(OrderedzSP_ncTX+min(OrderedzSP_ncTX(:))+1));
t_lim = [datenum(DateNum_IC_TimeAxis(STI_IC)) datenum(DateNum_IC_TimeAxis(ETI_IC))];
xlim(t_lim);
datetick('x','keeplimits');
ylabel('ncTX')
pause(0.1)


subplot(4,1,4)
OrderedzSP_sp = zSP_sp(:,MaxFiringRateOrder);
S = plotRasterplot(STI_IC,ETI_IC,fignum,DateNum_IC_TimeAxis,OrderedzSP_sp);
%S = plotBandplot(STI_IC,ETI_IC,k,DateNum_IC_TimeAxis,(OrderedzSP_sp+min(OrderedzSP_sp(:))+1));
t_lim = [datenum(DateNum_IC_TimeAxis(STI_IC)) datenum(DateNum_IC_TimeAxis(ETI_IC))];
xlim(t_lim);
datetick('x','keeplimits');
ylabel('SpikePower')
pause(0.1)
end

end

%% %% Now identify the order for 2x speed
i = 1;
j = 12; %(for 1x speed)

IndicesToLookAt = CCAllInterestingIndicesModified{i,j};

%Game playing:
%xLim1 = datetime(2020,11,24,14,50,0);
%xLim2 = datetime(2020,11,24,17,15,0);
%Overnight    
xLim1 = datetime(2020,11,24,23,0,0);
xLim2 = datetime(2020,11,25,9,0,0);
    
Indices_of_TwoXReplayEvents = TimeAxis(IndicesToLookAt);

a = find(Indices_of_TwoXReplayEvents>xLim1 & Indices_of_TwoXReplayEvents<xLim2);

Indices_of_TwoXReplayEvents = IndicesToLookAt(a);

TwoXRePlaySpikeOrder = zeros(length(Indices_of_TwoXReplayEvents),96);


sigmaF = .15;   
oSP = zSP_ncTX;
%oSP = zSP_sp;
    tic
    for j = 1:96
        oSP_j = oSP(:,j);
        smoothed_oSP_j = conv(oSP_j,exp(-(((-10:0.1:10).^2)./(2*sigmaF.^2))),'same');
        oSP(:,j) = smoothed_oSP_j;
    end
    toc
    
    
for k = 1:length(Indices_of_TwoXReplayEvents)
%for k = 50:70
STI_IC = Indices_of_TwoXReplayEvents(k) - 50*1; 
ETI_IC = Indices_of_TwoXReplayEvents(k) + 50*5;

%Find the order of firing based on time of max firing rate (try spikePower
%and ncTX separately and see which works better) within the X seconds of
%STCE
    
    TimeWindowToLookForPeak = 2;
    
    oSPi = oSP((Indices_of_TwoXReplayEvents(k)):(Indices_of_TwoXReplayEvents(k)+50*TimeWindowToLookForPeak),:);
    
    clear maxSPindex
    for j = 1:96
       oSP_j = oSPi(:,j);
       [a b] = max(oSP_j);
       b=b(1);
       maxSPindex(j) = b;
    end

    [a MaxFiringRateOrder] = sort(maxSPindex,'descend');
    
    [a FiringRateOrder] = sort(maxSPindex,'ascend');
    TwoXRePlaySpikeOrder(k,:) = FiringRateOrder;
    
%now plot the spiking data in order:
fignum = k*0+1;
figure(fignum)
set(fignum,'Position',[488 -134.2000 560 896.2000])
plotplot = 0;

if plotplot == 1
subplot(4,1,1)
OrderedThresholded_zSP = Thresholded_zSP(:,MaxFiringRateOrder);
S = plotRasterplot(STI_IC,ETI_IC,fignum,DateNum_IC_TimeAxis,OrderedThresholded_zSP);
%plot(DateNum_IC_TimeAxis(STI_IC:ETI_IC),mean(SP(STI_IC:ETI_IC,:),2))
t_lim = [datenum(DateNum_IC_TimeAxis(STI_IC)) datenum(DateNum_IC_TimeAxis(ETI_IC))];
xlim(t_lim);
datetick('x','keeplimits');
%ylabel('Mean Firing Rate')
pause(0.1)

subplot(4,1,2)
hold on
plot(DateNum_IC_TimeAxis(STI_IC:ETI_IC),Kalman_Data(1:2,STI_IC:ETI_IC))
plot([DateNum_IC_TimeAxis(Indices_of_TwoXReplayEvents(k)) DateNum_IC_TimeAxis(Indices_of_TwoXReplayEvents(k))],[-3 3],'g--')
axis tight; box off
datetick('x','keeplimits');
title('Kalman Filter Output')

subplot(4,1,3)
OrderedzSP_ncTX = zSP_ncTX(:,MaxFiringRateOrder);
S = plotRasterplot(STI_IC,ETI_IC,fignum,DateNum_IC_TimeAxis,OrderedzSP_ncTX);
%S = plotBandplot(STI_IC,ETI_IC,k,DateNum_IC_TimeAxis,(OrderedzSP_ncTX+min(OrderedzSP_ncTX(:))+1));
t_lim = [datenum(DateNum_IC_TimeAxis(STI_IC)) datenum(DateNum_IC_TimeAxis(ETI_IC))];
xlim(t_lim);
datetick('x','keeplimits');
ylabel('ncTX')
pause(0.1)


subplot(4,1,4)
OrderedzSP_sp = zSP_sp(:,MaxFiringRateOrder);
S = plotRasterplot(STI_IC,ETI_IC,fignum,DateNum_IC_TimeAxis,OrderedzSP_sp);
%S = plotBandplot(STI_IC,ETI_IC,k,DateNum_IC_TimeAxis,(OrderedzSP_sp+min(OrderedzSP_sp(:))+1));
t_lim = [datenum(DateNum_IC_TimeAxis(STI_IC)) datenum(DateNum_IC_TimeAxis(ETI_IC))];
xlim(t_lim);
datetick('x','keeplimits');
ylabel('SpikePower')
pause(0.1)
end

end


%% Calculate matching index between each replay event

%M is the total number of channels (96 here) in a given replay event (e.g.
%all of them).  There are a total of M(M-1)/2 possible cell pairs
% if "m" is the number of pairs that have the same order between a template
% and a frame, and "n" is the number of pairs with opposite order, then

% I = (m - n)/(m + n) 

% Where I is a matching index bounded by [-1 to 1], where 1 is a perfect
% match and -1 is an inverted match

%%%
%ActualGamePlaySpikeOrder are the 1x events.  Numbers 1 -- 89 are during
%the task performance; numbers 103 -- 187 are during overnight 

[numEvents numChannels] = size(ActualGamePlaySpikeOrder)
matchingmatrix = zeros(numEvents);
matchingmatrix1x = zeros(89,85);
tic
%for i = 1:numEvents
for i = 1:89
   template =  ActualGamePlaySpikeOrder(i,:);
   %for j = 1:numEvents
   for j = 1:85
       frame_to_match = ActualGamePlaySpikeOrder(j+102,:);
              
       Imatch = calculatematchingindex(template,frame_to_match);  
       matchingmatrix1x(i,j) = Imatch;
   end
   i
end
toc


% 
% % Generate a random control matching template
% clear controlImatch
% tic
% for i = 1:100000
%    template = randperm(96);
%    frame_to_match = randperm(96);
%    
%    Imatch = calculatematchingindex(template,frame_to_match);
%    controlImatch(i) = Imatch;   
% end
% toc
%%
clear differentFromControl
for i = 1:89
   thisDistribution = matchingmatrix1x(i,:);
   
   [h p] = ttest2(thisDistribution, controlImatch);
   differentFromControl(i) = p;
    
end
    
%save('MatchingVariables_Controls.mat','matchingmatrix','controlImatch')

%% Calculate matching index between daytime performance and 2x replay events


[numEvents_taskperformance numChannels] = size(ActualGamePlaySpikeOrder)
[numEvents_2xreplay numChannels] = size(TwoXRePlaySpikeOrder)

matchingmatrix2x = zeros(89,numEvents_2xreplay);

tic
%for i = 1:numEvents_taskperformance
for i = 1:89
   template =  ActualGamePlaySpikeOrder(i,:);
   for j = 1:numEvents_2xreplay
       frame_to_match = TwoXRePlaySpikeOrder(j,:);
              
       Imatch = calculatematchingindex(template,frame_to_match);  
       matchingmatrix2x(i,j) = Imatch;
   end
   i
end
toc


%% Play Overnight Activity ordered by a daytime task playing sequence
%First Load the 40-125 Hz LFP band for plotting

GammaBand = load('Band4_Array_2_SecondSession.mat');
GammaBand = GammaBand.Band4_array2;

%%
for zzz = 31 

%Game playing:
xLim1 = datetime(2020,11,24,14,50,0);
%xLim2 = datetime(2020,11,24,17,15,0);
%Overnight    
%xLim1 = datetime(2020,11,24,23,0,0);
xLim2 = datetime(2020,11,25,9,0,0);

UseThisOnesOrder = zzz;

Indices_of_OvernightEvents = TimeAxis(IndicesToLookAt);

a = find(Indices_of_OvernightEvents>xLim1 & Indices_of_OvernightEvents<xLim2);

Indices_of_OvernightEvents = IndicesToLookAt(a);

for k = UseThisOnesOrder

STI_IC = Indices_of_actualGamePlayEvents(k) - 50*4; 
ETI_IC = Indices_of_actualGamePlayEvents(k) + 50*8;

%Find the order of firing based on time of max firing rate (try spikePower
%and ncTX separately and see which works better) within the X seconds of
%STCE
    
    TimeWindowToLookForPeak = 3;
    
    oSPi = oSP((Indices_of_actualGamePlayEvents(k)+25):(Indices_of_actualGamePlayEvents(k)+50*TimeWindowToLookForPeak),:);
    
    clear maxSPindex
    for j = 1:96
       oSP_j = oSPi(:,j);
       [a b] = max(oSP_j);
       b=b(1);
       maxSPindex(j) = b;
    end

    [a MaxFiringRateOrder] = sort(maxSPindex,'descend');
    
    [a FiringRateOrder] = sort(maxSPindex,'ascend');
    ActualGamePlaySpikeOrder(k,:) = FiringRateOrder;
    
    Ordered_ozSP = oSP(:,MaxFiringRateOrder);
    oSPi = Ordered_ozSP((Indices_of_actualGamePlayEvents(k)+25):(Indices_of_actualGamePlayEvents(k)+50*TimeWindowToLookForPeak),:);
    
    clear maxSPindex
    for j = 1:96
       oSP_j = oSPi(:,j);
       [a b] = max(oSP_j);
       b=b(1);
       maxSPindex(j) = b;
    end
        
    p = polyfit(1:96,maxSPindex,1);
    pf = polyval(p,1:96);
    [h p] = corr((1:96)',maxSPindex','type','spearman');
    
%now plot the spiking data in order:
fignum = k;
figure(fignum)
clf
figure(fignum)
set(fignum,'Position',[50 -150 800 1000])
set(fignum,'DefaultAxesFontSize',16);

plotplot = 1;

if plotplot == 1

% Plot the EEG data for a minute around the ripple
Xlim1_eeg = TimeAxis(STI_IC) - minutes(1);
Xlim2_eeg = TimeAxis(ETI_IC) + minutes(1);

subplot(3,1,1)
h = plotEEGSpectrogram(Xlim1_eeg,Xlim2_eeg,fignum,EEGSpectrogramTimeAxis,EEG_Spectrogram_Frequency_Axis,EEG_Spectrogram);
title('EEG Spectrogram')
 
% subplot(4,1,2)
% h = plotBandplot(STI_IC,ETI_IC,fignum,DateNum_IC_TimeAxis,GammaBand,[2 4]);
% title('Gamma (40-125 Hz)')

subplot(3,1,2)
hold on
OrderedThresholded_zSP = Thresholded_zSP(:,MaxFiringRateOrder);
S = plotRasterplot(STI_IC,ETI_IC,fignum,DateNum_IC_TimeAxis,OrderedThresholded_zSP);
%plot(DateNum_IC_TimeAxis(STI_IC:ETI_IC),mean(SP(STI_IC:ETI_IC,:),2))
polyline = datenum(TimeAxis(Indices_of_actualGamePlayEvents(k)+25) + seconds(pf*0.02));
plot(polyline,1:96,'g-')
t_lim = [datenum(DateNum_IC_TimeAxis(STI_IC)) datenum(DateNum_IC_TimeAxis(ETI_IC))];
xlim(t_lim);
datetick('x','keeplimits');
%ylabel('Mean Firing Rate')
pause(0.1)

subplot(3,1,3)
hold on
plot(DateNum_IC_TimeAxis(STI_IC:ETI_IC),Kalman_Data(1:2,STI_IC:ETI_IC))
%plot([DateNum_IC_TimeAxis(Indices_of_actualGamePlayEvents(k)) DateNum_IC_TimeAxis(Indices_of_actualGamePlayEvents(k))],[-3 3],'g--')
axis tight; box off
%ylim([-2.5 2.5])
datetick('x','keeplimits');
title('Kalman Filter Output')
%
TimeOfEvent = Xlim1_eeg+minutes(1);
PrintFigTitle = ['Example Replay Spiking Order_template_' num2str(UseThisOnesOrder) '_ for event at ' num2str(datenum(TimeOfEvent)) '.svg'];
print(PrintFigTitle,'-dsvg','-painters')

end
    
end

% Plot some 1x examples:

%for k = 103:length(Indices_of_actualGamePlayEvents)
%for k = [106 117 139 157]

TheseIMatchess = (matchingmatrix1x(UseThisOnesOrder,:));
z = find(TheseIMatchess > 0.15);
z = z+102;
for  k = z
  
STI_IC = Indices_of_actualGamePlayEvents(k) - 50*0.5; 
ETI_IC = Indices_of_actualGamePlayEvents(k) + 50*4.5;


% Plot the EEG data for a minute around the ripple
Xlim1_eeg = TimeAxis(STI_IC) - minutes(1);
Xlim2_eeg = TimeAxis(ETI_IC) + minutes(1);


%now plot the spiking data in order:
fignum = 11;
figure(fignum)
clf
figure(fignum)
set(fignum,'Position',[50 -150 800 1000])
set(fignum,'DefaultAxesFontSize',16);

plotplot = 1;

if plotplot == 1

subplot(3,1,1)
h = plotEEGSpectrogram(Xlim1_eeg,Xlim2_eeg,fignum,EEGSpectrogramTimeAxis,EEG_Spectrogram_Frequency_Axis,EEG_Spectrogram);
title('EEG Spectrogram')

    
%%% Find the best fit line between plotted order and actual order:

%Find the order of firing based on time of max firing rate (try spikePower
%and ncTX separately and see which works better) within the X seconds of
%STCE
   
    TimeWindowToLookForPeak = 3;
    Ordered_ozSP = oSP(:,MaxFiringRateOrder);
    
    oSPi = Ordered_ozSP((Indices_of_actualGamePlayEvents(k)+25):(Indices_of_actualGamePlayEvents(k)+50*TimeWindowToLookForPeak),:);
    
    clear maxSPindex
    for j = 1:96
       oSP_j = oSPi(:,j);
       [a b] = max(oSP_j);
       b=b(1);
       maxSPindex(j) = b;
    end

    p = polyfit(1:96,maxSPindex,1);
    pf = polyval(p,1:96);
    [h p] = corr((1:96)',maxSPindex','type','spearman');
    
subplot(3,1,2)
hold on
%OrderedThresholded_zSP = Thresholded_zSP(:,MaxFiringRateOrder);
OrderedThresholded_zSP = Thresholded_zSP(:,MaxFiringRateOrder);
S = plotRasterplot(STI_IC,ETI_IC,fignum,DateNum_IC_TimeAxis,OrderedThresholded_zSP);
%plot(DateNum_IC_TimeAxis(STI_IC:ETI_IC),mean(SP(STI_IC:ETI_IC,:),2))
polyline = datenum(TimeAxis(Indices_of_actualGamePlayEvents(k)+25) + seconds(pf*0.02));
plot(polyline,1:96,'g-')
t_lim = [datenum(DateNum_IC_TimeAxis(STI_IC)) datenum(DateNum_IC_TimeAxis(ETI_IC))];
xlim(t_lim);
datetick('x','keeplimits');
%ylabel('Mean Firing Rate')
pause(0.1)

subplot(3,1,3)
hold on
plot(DateNum_IC_TimeAxis(STI_IC:ETI_IC),Kalman_Data(1:2,STI_IC:ETI_IC))
%plot([DateNum_IC_TimeAxis(Indices_of_actualGamePlayEvents(k)) DateNum_IC_TimeAxis(Indices_of_actualGamePlayEvents(k))],[-3 3],'g--')
text(DateNum_IC_TimeAxis(Indices_of_actualGamePlayEvents(k)),1,['I = ' num2str(matchingmatrix1x(UseThisOnesOrder,k-102))])
text(DateNum_IC_TimeAxis(Indices_of_actualGamePlayEvents(k)),0.5,['p = ' num2str(p)])
axis tight; box off
t_lim = [datenum(DateNum_IC_TimeAxis(STI_IC)) datenum(DateNum_IC_TimeAxis(ETI_IC))];
xlim(t_lim);
%ylim([-2.5 2.5])
datetick('x','keeplimits');
title('Kalman Filter Output')
pause(1)

% 
end

TimeOfEvent = Xlim1_eeg+minutes(1);
PrintFigTitle = ['Example Replay Spiking Order_template_' num2str(UseThisOnesOrder) '_ at 1x event at ' num2str(datenum(TimeOfEvent)) '.svg'];
print(PrintFigTitle,'-dsvg','-painters')

end


% Plot some 2x examples:
%for k = 1:length(Indices_of_TwoXReplayEvents)

TheseIMatchess = (matchingmatrix2x(UseThisOnesOrder,:));
z = find(TheseIMatchess > 0.15);

for k = z
    
  
STI_IC = Indices_of_TwoXReplayEvents(k) - 50*0.5; 
ETI_IC = Indices_of_TwoXReplayEvents(k) + 50*2.5;


% Plot the EEG data for a minute around the ripple
Xlim1_eeg = TimeAxis(STI_IC) - minutes(1);
Xlim2_eeg = TimeAxis(ETI_IC) + minutes(1);


%now plot the spiking data in order:
fignum = 22;
figure(fignum)
clf
figure(fignum)
set(fignum,'Position',[50 -150 800 1000])
set(fignum,'DefaultAxesFontSize',16);

plotplot = 1;

if plotplot == 1

subplot(3,1,1)
h = plotEEGSpectrogram(Xlim1_eeg,Xlim2_eeg,fignum,EEGSpectrogramTimeAxis,EEG_Spectrogram_Frequency_Axis,EEG_Spectrogram);
title('EEG Spectrogram')

   
%%% Find the best fit line between plotted order and actual order:

%Find the order of firing based on time of max firing rate (try spikePower
%and ncTX separately and see which works better) within the X seconds of
%STCE
   
    TimeWindowToLookForPeak = 2;
    Ordered_ozSP = oSP(:,MaxFiringRateOrder);
    
    oSPi = Ordered_ozSP((Indices_of_TwoXReplayEvents(k)+0):(Indices_of_TwoXReplayEvents(k)+50*TimeWindowToLookForPeak),:);
    
    clear maxSPindex
    for j = 1:96
       oSP_j = oSPi(:,j);
       [a b] = max(oSP_j);
       b=b(1);
       maxSPindex(j) = b;
    end

    p = polyfit(1:96,maxSPindex,1);
    pf = polyval(p,1:96);
    [h p] = corr((1:96)',maxSPindex','type','spearman');
    
    
subplot(3,1,2)
hold on
%OrderedThresholded_zSP = Thresholded_zSP(:,MaxFiringRateOrder);
OrderedThresholded_zSP = Thresholded_zSP(:,MaxFiringRateOrder);
S = plotRasterplot(STI_IC,ETI_IC,fignum,DateNum_IC_TimeAxis,OrderedThresholded_zSP);
%plot(DateNum_IC_TimeAxis(STI_IC:ETI_IC),mean(SP(STI_IC:ETI_IC,:),2))
polyline = datenum(TimeAxis(Indices_of_TwoXReplayEvents(k)+0) + seconds(pf*0.02));
plot(polyline,1:96,'g-')
t_lim = [datenum(DateNum_IC_TimeAxis(STI_IC)) datenum(DateNum_IC_TimeAxis(ETI_IC))];
xlim(t_lim);
datetick('x','keeplimits');
%ylabel('Mean Firing Rate')
pause(0.1)

subplot(3,1,3)
hold on
plot(DateNum_IC_TimeAxis(STI_IC:ETI_IC),Kalman_Data(1:2,STI_IC:ETI_IC))
%plot([DateNum_IC_TimeAxis(Indices_of_TwoXReplayEvents(k)) DateNum_IC_TimeAxis(Indices_of_TwoXReplayEvents(k))],[-3 3],'g--')
text(DateNum_IC_TimeAxis(Indices_of_TwoXReplayEvents(k)),1,['I = ' num2str(matchingmatrix2x(UseThisOnesOrder,k))])
text(DateNum_IC_TimeAxis(Indices_of_TwoXReplayEvents(k)),0.5,['p = ' num2str(p)])

axis tight; box off
t_lim = [datenum(DateNum_IC_TimeAxis(STI_IC)) datenum(DateNum_IC_TimeAxis(ETI_IC))];
xlim(t_lim);
%ylim([-2.5 2.5])
datetick('x','keeplimits');
title('Kalman Filter Output')
pause(1)

% 


TimeOfEvent = Xlim1_eeg+minutes(1);
PrintFigTitle = ['S2 Example Replay Spiking Order_template_' num2str(UseThisOnesOrder) '_ at 2x event at ' num2str(datenum(TimeOfEvent)) '.svg'];
print(PrintFigTitle,'-dsvg','-painters')

%Title_Details = sprintf(', %s to %s',datetime(Xlim1,'Format','HH-mm-ss'),datetime(Xlim2,'Format','HH-mm-ss'));
%PrintFigTitle = ['An Example of Ripple at ' Title_Details ' play at ' num2str(templateCompressionFactor(k)) 'x.svg'];

end

end




end