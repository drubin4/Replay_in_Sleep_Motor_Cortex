
%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021

%% Evaluate task performance the day after sleeping

%% Open intracranial neural data as well as surface EEG data:

clear all
close all
clc

%% Add all directories to path to allow for ease of access to files

a = pwd;

if strcmp(a(end-3:end),'ipts')
    cd ..
    uponedir = pwd;
    addpath(genpath(uponedir));
end
if strcmp(a(end-3:end),'sion')
    uponedir = pwd;
    addpath(genpath(uponedir));
end
    
%% First open the intracranial data set of interest:

%DataSetToOpen = 'ncTX_Array_1.mat';
DataSetToOpen = 'ncTX_Array_2.mat';
%DataSetToOpen = 'SpikePower_Array_1.mat';
%DataSetToOpen = 'SpikePower_Array_2.mat';

load(DataSetToOpen)

if strcmp(DataSetToOpen,'ncTX_Array_1.mat')
SP = single(ncTX_array1);                   %Array 1 is lateral array
NeuralCodeSource = 'ncTX_array1';
clear ncTX_array1
end
if strcmp(DataSetToOpen,'ncTX_Array_2.mat')
SP = single(ncTX_array2);                  %Array 2 is medial array
NeuralCodeSource = 'ncTX_array2';
clear ncTX_array2
end
if strcmp(DataSetToOpen,'SpikePower_Array_1.mat')
SP = SpikePower_array1;
NeuralCodeSource = 'SpikePower_array1';
clear SpikePower_array1
end
if strcmp(DataSetToOpen,'SpikePower_Array_2.mat')
SP = SpikePower_array2;
NeuralCodeSource = 'SpikePower_array2';
clear SpikePower_array2
end

%% Next open the EEG File:

load('Synced_EEG_Data.mat')

%Time Axis for EEG is in 512 Hz
%Time Axis for EEG spectrogram is in Seconds:
EEGSpectrogramTimeAxis = datenum(EEG_Spectrogram_DateTimeAxis);


T_start_EEG = datetime(2020,11,23,18,0,0);
T_end_EEG = datetime(2020,11,25,11,30,0);


% Plot The EEG Spectrogram (just cuz, for display purposes)
fignum = 3;
h = plotEEGSpectrogram(T_start_EEG,T_end_EEG,fignum,EEGSpectrogramTimeAxis,EEG_Spectrogram_Frequency_Axis,EEG_Spectrogram);

% Z-score and threshold the IC data for plotting Rasters (these are also for display; no analysis is done on these data):
clear zSP

%z-score each channel
for i = 1:96
    
   SPi = SP(:,i);  
   zSP(:,i) = (SPi - mean(SPi))/std(SPi);   
    
end

%Threshold at > 1 SD above mean
Thresholded_zSP = zSP > 1;
Thresholded_zSP = double(Thresholded_zSP);

%Pad with zeros at discontinuities for plotting to prevent weird streaks on the Raster plot:
TimeAxisDiff = diff(TimeAxis);
jumps = find(TimeAxisDiff>seconds(1));

Thresholded_zSP(jumps,:) = 0;
Thresholded_zSP(jumps+1,:) = 0;

% Pull in all of the features to be used by the KalmanModel

load('KalmanModelData.mat')
A_KM = KalmanModel.A;
H_KM = KalmanModel.H;
K_KM = KalmanModel.K;

ncTXInds = KalmanModel.ncTXInds;
spikePowerInds = KalmanModel.spikePowerInds;

load('ncTX_Both_Array.mat')
ncTX_KM = ncTX_Both_Arrays(:,ncTXInds);
clear ncTX_Both_Arrays

load('SpikePower_Both_Arrays.mat')
spikePower_KM = SpikePower_Both_Arrays(:,spikePowerInds);
clear SpikePower_Both_Arrays

KM_features = (cat(2,ncTX_KM,spikePower_KM));   %% KM_features is the duration of the session x 40 channels (13 nctx and 27 spike power) used by the Kalman filter model to generate the simualated movement


%%
%% Look at a Simon Block, and draw the trajectory of the cursor on each round of the game

DateNum_IC_TimeAxis = datenum(TimeAxis);

GamePerformance = zeros(4,3);

for WhichSimonBlock = [1 2 4]

    [GamePosition GameSequences GameClock GameTrialOutcomes GameNewTrialTimes GameNewTrialEndTimes Xlim1 Xlim2] = ExtractGameDataNextDay(WhichSimonBlock)
    
 StartTimeIndex_IC = find(TimeAxis > Xlim1);
 StartTimeIndex_IC = StartTimeIndex_IC(1);

EndTimeIndex_IC = find(TimeAxis > Xlim2);
EndTimeIndex_IC = EndTimeIndex_IC(1)-1;

DateNum_IC_TimeAxis = datenum(TimeAxis);

t_lim = [datenum(Xlim1) datenum(Xlim2)];
%%%%%Run the data through the Kalman Decoder Model
CursorPositionOutput = PlayNeuralDataThroughKalmanEML(KM_features,StartTimeIndex_IC,EndTimeIndex_IC,A_KM,H_KM,K_KM);

%%%Plot
fignum = WhichSimonBlock; 
figure(fignum)
clf
figure(fignum)
set(fignum,'Position',[0 -100 1600 800])

% 
% subplot(4,16,1:16)
% h = plotEEGSpectrogram(Xlim1,Xlim2,fignum,EEGSpectrogramTimeAxis,EEG_Spectrogram_Frequency_Axis,EEG_Spectrogram);

subplot(4,16,17:32)
S = plotRasterplot(StartTimeIndex_IC,EndTimeIndex_IC,fignum,DateNum_IC_TimeAxis,Thresholded_zSP);


subplot(4,16,33:48)
hold on

for j = 1:length(GameNewTrialTimes)
   pos = zeros(1,4);
   pos(1) = GameNewTrialTimes(j);
   pos(2) = -0.5;
   pos(3) = GameNewTrialEndTimes(j) - GameNewTrialTimes(j);
   pos(4) = 1;  
   if GameTrialOutcomes(j) == 1
       h = rectangle('Position',pos,'EdgeColor','none','FaceColor',[0.5 1 0.5]);
   else
       h = rectangle('Position',pos,'EdgeColor','none');
       h.FaceColor = [1 0.65 0.65];
   end   

end

plot(GameClock,GamePosition)

axis tight;
xlim(t_lim)
ylim([-0.4 0.4])
datetick('x','keeplimits');
ylabel('Cursor Position');
title('Cursor Position')




subplot(4,16,49:64)

hold on
for j = 1:length(GameNewTrialTimes)
   pos = zeros(1,4);
   pos(1) = GameNewTrialTimes(j);
   pos(2) = -3;
   pos(3) = GameNewTrialEndTimes(j) - GameNewTrialTimes(j);
   pos(4) = 6;  
   if sum((GameSequences(j,:) == [2 4 1 3])) == 4
       h = rectangle('Position',pos,'EdgeColor','none','FaceColor',[0.85 0.85 0.85]);
   else
       h = rectangle('Position',pos,'EdgeColor','none');
       h.FaceColor = [0.65 0.65 0.65];
   end
end

plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),CursorPositionOutput(1:2,:))

axis tight;
xlim(t_lim)
%ylim([-220 220])
datetick('x','keeplimits');
ylabel('Magnitude');
title('Raw Kalman Model Output')



Title_Details = sprintf(', %s to %s',datetime(Xlim1,'Format','HH-mm-ss'),datetime(Xlim2,'Format','HH-mm-ss'));
suptitle(['EEG, MEA, Cursor Position, and Kalman Filter Data' Title_Details])

PrintFigTitle = ['EEG, MEA, Cursor Position, and Kalman Filter Data' Title_Details  ' ' NeuralCodeSource '.pdf'];



fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 18 11];
fig.PaperSize = [18 11];
if printfigs == 1
print(PrintFigTitle,'-bestfit','-dpdf')
end

%Calculate performance statistics on game:

for j = 1:length(GameNewTrialTimes)
  
   if GameTrialOutcomes(j) == 1
        GamePerformance(WhichSimonBlock,1) = GamePerformance(WhichSimonBlock,1)+1;
        if sum((GameSequences(j,:) == [2 4 1 3])) == 4
            GamePerformance(WhichSimonBlock,2) = GamePerformance(WhichSimonBlock,2)+1;
        else
            GamePerformance(WhichSimonBlock,3) = GamePerformance(WhichSimonBlock,3)+1;
        end
   end   
   
end


end

%% Graph game performance:
GPerm = GamePerformance;
GPerm(3,:) = [];

OverallPerformance = GPerm(:,1)./16;
TargetPerformance = GPerm(:,2)./12;
DistractorPerformance = GPerm(:,3)./4;

CumulativeOverallPerformance = sum(GPerm(:,1))/(16*3);
CumulativeTargetPerformance = sum(GPerm(:,2))/(12*3);
CumulativeDistrctorPerformance = sum(GPerm(:,3))/(4*3);

NormalizedPerformance = [OverallPerformance TargetPerformance DistractorPerformance];
[h pTtest] = ttest2(TargetPerformance,DistractorPerformance);

pAnova1 = anova1(NormalizedPerformance);

figure(222)
clf
figure(222)
hold on

plot(1:3,TargetPerformance,'gs-')
plot(1:3,DistractorPerformance,'rs-')
plot(1:3,OverallPerformance,'bs-')

legend('Target','Distractor','Overall','Location','SouthWest')

ylim([0 1])
ylabel('Performance')
xlabel('Block #')

text(1.75,0.25,['Overall Performance: ' num2str(CumulativeOverallPerformance)])
text(1.75,0.20,['Target Performance: ' num2str(CumulativeTargetPerformance)])
text(1.75,0.15,['Distractor Performance: ' num2str(CumulativeDistrctorPerformance)])
text(1.75,0.1,['ANOVA p-value: ' num2str(pAnova1)])
text(1.75,0.05,['t test Target vs. Distractor p-value: ' num2str(pTtest)])

%print('After_Sleeping_Overnight_TaskPerformance.pdf','-bestfit','-dpdf')

%% Generate performance stats from first day:

%% Look at a single Simon Block, and draw the trajectory of the cursor on each round of the game

GamePerformanceBeforeSleep = zeros(10,3);

for WhichSimonBlock = 1:10     

 
    [GamePosition GameSequences GameClock GameTrialOutcomes GameNewTrialTimes GameNewTrialEndTimes Xlim1 Xlim2] = ExtractGameData(WhichSimonBlock)
    
 StartTimeIndex_IC = find(TimeAxis > Xlim1);
 StartTimeIndex_IC = StartTimeIndex_IC(1);

EndTimeIndex_IC = find(TimeAxis > Xlim2);
EndTimeIndex_IC = EndTimeIndex_IC(1)-1;

DateNum_IC_TimeAxis = datenum(TimeAxis);

t_lim = [datenum(Xlim1) datenum(Xlim2)];
%%%%%Run the data through the Kalman Decoder Model
CursorPositionOutput = PlayNeuralDataThroughKalmanEML(KM_features,StartTimeIndex_IC,EndTimeIndex_IC,A_KM,H_KM,K_KM);

%%%Plot
fignum = 200103 + WhichSimonBlock; 
figure(fignum)
clf
figure(fignum)
set(fignum,'Position',[0 -100 1600 800])

subplot(4,16,33:48)
hold on

for j = 1:length(GameNewTrialTimes)
   pos = zeros(1,4);
   pos(1) = GameNewTrialTimes(j);
   pos(2) = -0.5;
   pos(3) = GameNewTrialEndTimes(j) - GameNewTrialTimes(j);
   pos(4) = 1;  
   if GameTrialOutcomes(j) == 1
       h = rectangle('Position',pos,'EdgeColor','none','FaceColor',[0.5 1 0.5]);
   else
       h = rectangle('Position',pos,'EdgeColor','none');
       h.FaceColor = [1 0.65 0.65];
   end   

end

plot(GameClock,GamePosition)

axis tight;
xlim(t_lim)
ylim([-0.4 0.4])
datetick('x','keeplimits');
ylabel('Cursor Position');
title('Cursor Position')




subplot(4,16,49:64)

hold on
for j = 1:length(GameNewTrialTimes)
   pos = zeros(1,4);
   pos(1) = GameNewTrialTimes(j);
   pos(2) = -3;
   pos(3) = GameNewTrialEndTimes(j) - GameNewTrialTimes(j);
   pos(4) = 6;  
   if sum((GameSequences(j,:) == [2 4 1 3])) == 4
       h = rectangle('Position',pos,'EdgeColor','none','FaceColor',[0.85 0.85 0.85]);
   else
       h = rectangle('Position',pos,'EdgeColor','none');
       h.FaceColor = [0.65 0.65 0.65];
   end
end

plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),CursorPositionOutput(1:2,:))

axis tight;
xlim(t_lim)
%ylim([-220 220])
datetick('x','keeplimits');
ylabel('Magnitude');
title('Raw Kalman Model Output')

%Calculate performance statistics on game:

for j = 1:length(GameNewTrialTimes)
  
   if GameTrialOutcomes(j) == 1
        GamePerformanceBeforeSleep(WhichSimonBlock,1) = GamePerformanceBeforeSleep(WhichSimonBlock,1)+1;
        if sum((GameSequences(j,:) == [2 4 1 3])) == 4
            GamePerformanceBeforeSleep(WhichSimonBlock,2) = GamePerformanceBeforeSleep(WhichSimonBlock,2)+1;
        else
            GamePerformanceBeforeSleep(WhichSimonBlock,3) = GamePerformanceBeforeSleep(WhichSimonBlock,3)+1;
        end
   end   
   
end




end

%% Graph game performance:

OverallPerformance = GamePerformanceBeforeSleep(:,1)./16;
TargetPerformance = GamePerformanceBeforeSleep(:,2)./12;
DistractorPerformance = GamePerformanceBeforeSleep(:,3)./4;

CumulativeOverallPerformance = sum(GamePerformanceBeforeSleep(:,1))/(16*10);
CumulativeTargetPerformance = sum(GamePerformanceBeforeSleep(:,2))/(12*10);
CumulativeDistrctorPerformance = sum(GamePerformanceBeforeSleep(:,3))/(4*10);

NormalizedPerformance = [OverallPerformance TargetPerformance DistractorPerformance];
[h pTtest] = ttest2(TargetPerformance,DistractorPerformance);

pAnova1 = anova1(NormalizedPerformance);

figure(224)
clf
figure(224)
hold on

plot(1:10,TargetPerformance,'gs-')
plot(1:10,DistractorPerformance,'rs-')
plot(1:10,OverallPerformance,'bs-')

legend('Target','Distractor','Overall','Location','NorthWest')

ylim([0 1])
ylabel('Performance')
xlabel('Block #')

text(2,0.25,['Overall Performance: ' num2str(CumulativeOverallPerformance)])
text(2,0.20,['Target Performance: ' num2str(CumulativeTargetPerformance)])
text(2,0.15,['Distractor Performance: ' num2str(CumulativeDistrctorPerformance)])
text(2,0.1,['ANOVA p-value: ' num2str(pAnova1)])
text(2,0.05,['t test Target vs. Distractor p-value: ' num2str(pTtest)])

%print('Day2_TaskPerformance.pdf','-bestfit','-dpdf')

%% Graph next day game performance:
GPerm = GamePerformance;
GPerm(3,:) = [];

OverallPerformance = GPerm(:,1)./16;
TargetPerformance = GPerm(:,2)./12;
DistractorPerformance = GPerm(:,3)./4;

CumulativeOverallPerformance = sum(GPerm(:,1))/(16*3);
CumulativeTargetPerformance = sum(GPerm(:,2))/(12*3);
CumulativeDistrctorPerformance = sum(GPerm(:,3))/(4*3);

NormalizedPerformance = [OverallPerformance TargetPerformance DistractorPerformance];
[h pTtest] = ttest2(TargetPerformance,DistractorPerformance);

pAnova1 = anova1(NormalizedPerformance);

figure(222)
clf
figure(222)
hold on

plot(1:3,TargetPerformance,'gs-')
plot(1:3,DistractorPerformance,'rs-')
plot(1:3,OverallPerformance,'bs-')

legend('Target','Distractor','Overall','Location','SouthWest')

ylim([0 1])
ylabel('Performance')
xlabel('Block #')

text(1.75,0.25,['Overall Performance: ' num2str(CumulativeOverallPerformance)])
text(1.75,0.20,['Target Performance: ' num2str(CumulativeTargetPerformance)])
text(1.75,0.15,['Distractor Performance: ' num2str(CumulativeDistrctorPerformance)])
text(1.75,0.1,['ANOVA p-value: ' num2str(pAnova1)])
text(1.75,0.05,['t test Target vs. Distractor p-value: ' num2str(pTtest)])

%print('After_Sleeping_Overnight_TaskPerformance.pdf','-bestfit','-dpdf')


%% Graph both days performance


OverallPerformanceD1 = GamePerformanceBeforeSleep(:,1)./16;
TargetPerformanceD1 = GamePerformanceBeforeSleep(:,2)./12;
DistractorPerformanceD1 = GamePerformanceBeforeSleep(:,3)./4;

CumulativeOverallPerformanceD1 = sum(GamePerformanceBeforeSleep(:,1))/(16*10);
CumulativeTargetPerformanceD1 = sum(GamePerformanceBeforeSleep(:,2))/(12*10);
CumulativeDistrctorPerformanceD1 = sum(GamePerformanceBeforeSleep(:,3))/(4*10);

NormalizedPerformanceD1 = [OverallPerformanceD1 TargetPerformanceD1 DistractorPerformanceD1];
[h pTtestD1] = ttest2(TargetPerformanceD1,DistractorPerformanceD1);

pAnova1D1 = anova1(NormalizedPerformanceD1);


OverallPerformanceD2 = GPerm(:,1)./16;
TargetPerformanceD2 = GPerm(:,2)./12;
DistractorPerformanceD2 = GPerm(:,3)./4;

CumulativeOverallPerformanceD2 = sum(GPerm(:,1))/(16*3);
CumulativeTargetPerformanceD2 = sum(GPerm(:,2))/(12*3);
CumulativeDistrctorPerformanceD2 = sum(GPerm(:,3))/(4*3);

NormalizedPerformanceD2 = [OverallPerformanceD2 TargetPerformanceD2 DistractorPerformanceD2];
[h pTtestD2] = ttest2(TargetPerformanceD2,DistractorPerformanceD2);
pAnova1D2 = anova1(NormalizedPerformanceD2);



figure(225)
clf
figure(225)
set(225,'Position',[173   212   1200   420])

subplot(1,2,1)
hold on


plot(1:10,TargetPerformanceD1,'gs-')
plot(1:10,DistractorPerformanceD1,'rs-')
plot(1:10,OverallPerformanceD1,'bs-')

legend('Target','Distractor','Overall','Location','NorthWest')

ylim([0 1])
ylabel('Performance')
xlabel('Block #')

text(2,0.25,['Overall Performance: ' num2str(CumulativeOverallPerformanceD1)])
text(2,0.20,['Target Performance: ' num2str(CumulativeTargetPerformanceD1)])
text(2,0.15,['Distractor Performance: ' num2str(CumulativeDistrctorPerformanceD1)])
text(2,0.1,['ANOVA p-value: ' num2str(pAnova1D1)])
text(2,0.05,['t test Target vs. Distractor p-value: ' num2str(pTtestD1)])


subplot(1,2,2)
hold on

plot(1:3,TargetPerformanceD2,'gs-')
plot(1:3,DistractorPerformanceD2,'rs-')
plot(1:3,OverallPerformanceD2,'bs-')

legend('Target','Distractor','Overall','Location','SouthWest')

ylim([0 1])
ylabel('Performance')
xlabel('Block #')

text(1.65,0.25,['Overall Performance: ' num2str(CumulativeOverallPerformanceD2)])
text(1.65,0.20,['Target Performance: ' num2str(CumulativeTargetPerformanceD2)])
text(1.65,0.15,['Distractor Performance: ' num2str(CumulativeDistrctorPerformanceD2)])
text(1.65,0.1,['ANOVA p-value: ' num2str(pAnova1D2)])
text(1.65,0.05,['t test Target vs. Distractor p-value: ' num2str(pTtestD2)])

%Compare before and after sleep

[h pTtestOverall] = ttest2(OverallPerformanceD1,OverallPerformanceD2)
[h pTtestTarget] = ttest2(TargetPerformanceD1,TargetPerformanceD2)
[h pTtestDistractor] = ttest2(DistractorPerformanceD1,DistractorPerformanceD2)

suptitle({'Game Performance Before and After Sleep',['Overall p: ' num2str(pTtestOverall) '*, Target p: ' num2str(pTtestTarget) '*, Distractor p: ' num2str(pTtestDistractor)]})



fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 18 11];
fig.PaperSize = [18 11];
print('BothDays_TaskPerformance.pdf','-bestfit','-dpdf')
