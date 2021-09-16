%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021


%% The goal of this script will be to open neural data, reoncstruct the model used during the actual session, 
% and play the neural activity into the model to simulate the movement response.

%% Open intracranial neural data as well as surface EEG data:

clear all
close all
clc

printfigs = 0;       %Do (1) or do not (0) print figures as the analysis progresses
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

% First open an intracranial data set of interest (although technically this is just for plotting as
%the Kalman model uses data from both arrays and both spike power and NCTX)
%Array 1 is lateral array
%Array 2 is medial array

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
% Next open the EEG File: Also just for plotting here

load('Synced_EEG_Data.mat')

%Time Axis for EEG is in 512 Hz
RawEEGData;
EEG_Spectrogram_DateTimeAxis;

%Time Axis for EEG spectrogram is in Seconds:
EEG_Spectrogram;
EEG_Spectrogram_Frequency_Axis;
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

%Pad with zeros at discontinuities for plotting to prevent streaks on the Raster plot:
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

load('ncTX_Both_Arrays.mat')
ncTX_KM = ncTX_Both_Arrays(:,ncTXInds);
clear ncTX_Both_Arrays

load('SpikePower_Both_Arrays.mat')
spikePower_KM = SpikePower_Both_Arrays(:,spikePowerInds);
clear SpikePower_Both_Arrays

KM_features = (cat(2,ncTX_KM,spikePower_KM));   %% KM_features is the duration of the session x 40 channels (13 nctx and 27 spike power) used by the Kalman filter model to generate the simualated movement

%% Try running some data through the Kalman model to see what happens

%Define a brief segment of time
 Xlim1 = datetime(2020,11,24,16,57,30);
 Xlim2 = datetime(2020,11,24,17,0,30);

StartTimeIndex_IC = find(TimeAxis > Xlim1);
StartTimeIndex_IC = StartTimeIndex_IC(1);

EndTimeIndex_IC = find(TimeAxis > Xlim2);
EndTimeIndex_IC = EndTimeIndex_IC(1);

DateNum_IC_TimeAxis = datenum(TimeAxis);

%%%%%Run the neural data through the Kalman Decoder Model
%%%This code should really only use a few minutes of data at a time because
%%%the open-loop Kalman model is normalized in real time
CursorPositionOutput = PlayNeuralDataThroughKalmanEML(KM_features,StartTimeIndex_IC,EndTimeIndex_IC,A_KM,H_KM,K_KM);

%%Plot the data

figure(200103)
clf
figure(200103)
set(200103,'Position',[50 100 1000 600])

subplot(3,1,1)
h = plotEEGSpectrogram(Xlim1,Xlim2,200103,EEGSpectrogramTimeAxis,EEG_Spectrogram_Frequency_Axis,EEG_Spectrogram);

subplot(3,1,2)
S = plotRasterplot(StartTimeIndex_IC,EndTimeIndex_IC,200103,DateNum_IC_TimeAxis,Thresholded_zSP);

subplot(3,1,3)

hold on
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),CursorPositionOutput(1:2,:))

axis tight;
datetick('x','keeplimits');
ylabel('Magnitude');
title('Raw Kalman Model Output')

Title_Details = sprintf(', %s to %s',datetime(Xlim1,'Format','HH-mm-ss'),datetime(Xlim2,'Format','HH-mm-ss'));
suptitle(['EEG, MEA, and Predicted Cursor Position' Title_Details])

PrintFigTitle = ['EEG, MEA, and Predicted Cursor Position' Title_Details ' ' NeuralCodeSource '.pdf'];

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 18 11];
fig.PaperSize = [18 11];

if printfigs == 1
%print(PrintFigTitle,'-bestfit','-dpdf')
end

%% Identify which time segments correspond to different activity periods:


T_start_EEG = datetime(2020,11,24,0,0,0);
T_end_EEG = datetime(2020,11,25,9,0,0);


StartTimeIndex_IC = find(TimeAxis > T_start_EEG);
StartTimeIndex_IC = StartTimeIndex_IC(1);
EndTimeIndex_IC = find(TimeAxis > T_end_EEG);
EndTimeIndex_IC = EndTimeIndex_IC(1);
DateNum_IC_TimeAxis = datenum(TimeAxis);

%% Cycle through each Block of the matching game, and record the the neural activity and timing of each trial
% taking special note of which of the trials were successfull vs. not and which were target vs. distractor


for WhichSimonBlock = 1:10      

    [GamePosition GameSequences GameClock GameTrialOutcomes GameNewTrialTimes GameNewTrialEndTimes Xlim1 Xlim2] = ExtractGameData(WhichSimonBlock)
    
 StartTimeIndex_IC = find(TimeAxis > Xlim1);
 StartTimeIndex_IC = StartTimeIndex_IC(1);

EndTimeIndex_IC = find(TimeAxis > Xlim2);
EndTimeIndex_IC = EndTimeIndex_IC(1)-1;

DateNum_IC_TimeAxis = datenum(TimeAxis);

%%%%%Run the data through the Kalman Decoder Model
CursorPositionOutput = PlayNeuralDataThroughKalmanEML(KM_features,StartTimeIndex_IC,EndTimeIndex_IC,A_KM,H_KM,K_KM);


%%%Plot
figure(200103+WhichSimonBlock)
clf
figure(200103+WhichSimonBlock)
set(200103+WhichSimonBlock,'Position',[50 -300 1000 800])
t_lim = [datenum(Xlim1) datenum(Xlim2)];

subplot(4,1,1)
h = plotEEGSpectrogram(Xlim1,Xlim2,200103+WhichSimonBlock,EEGSpectrogramTimeAxis,EEG_Spectrogram_Frequency_Axis,EEG_Spectrogram);

subplot(4,1,2)
S = plotRasterplot(StartTimeIndex_IC,EndTimeIndex_IC,200103+WhichSimonBlock,DateNum_IC_TimeAxis,Thresholded_zSP);

subplot(4,1,3)

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
xlabel('Time');
title('Cursor Position')


subplot(4,1,4)

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
xlabel('Time');
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


end


%%
%% Look at a single Simon Block, and draw the trajectory of the cursor on each round of the game

GamePerformance = zeros(10,3);

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


subplot(4,16,1:16)
h = plotEEGSpectrogram(Xlim1,Xlim2,fignum,EEGSpectrogramTimeAxis,EEG_Spectrogram_Frequency_Axis,EEG_Spectrogram);

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

OverallPerformance = GamePerformance(:,1)./16;
TargetPerformance = GamePerformance(:,2)./12;
DistractorPerformance = GamePerformance(:,3)./4;

CumulativeOverallPerformance = sum(GamePerformance(:,1))/(16*10);
CumulativeTargetPerformance = sum(GamePerformance(:,2))/(12*10);
CumulativeDistrctorPerformance = sum(GamePerformance(:,3))/(4*10);

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

print('TaskPerformance.pdf','-bestfit','-dpdf')


%%


 %% Now push all the data from throughout the full session into the Kallman model, searching for evidence of learning of the target patter in the neural activity
 % Break the data into ~1 minute blocks, normalize the data, feed it into
 % Kalman model.
 
 % Then, with the Kalman model complete throughout, need to
 % do pattern matching to find where the target pattern is demonstrated -- more on that later.
 
 timeStart = datetime(2020,11,23,23,0,0);
 timeEnd = datetime(2020,11,25,9,0,0);
 
 duration_in_minutes = minutes(timeEnd - timeStart);
 analysis_interval = 1;                                  %unit is minutes
 number_of_analysis_segments = ceil(duration_in_minutes/analysis_interval);
 
DateNum_IC_TimeAxis = datenum(TimeAxis);

Kalman_Data = zeros(3,length(DateNum_IC_TimeAxis));

for i = 1:number_of_analysis_segments
     
     segmentTimeStart = timeStart + minutes(analysis_interval)*(i-1);
     segmentTimeEnd = segmentTimeStart+minutes(2);
    
     %%% First find if there is IC data present during the time interval
     
    StartTimeIndex_IC = find(TimeAxis > segmentTimeStart);
    StartTimeIndex_IC = StartTimeIndex_IC(1);

    EndTimeIndex_IC = find(TimeAxis > segmentTimeEnd);
    EndTimeIndex_IC = EndTimeIndex_IC(1);
    
    if EndTimeIndex_IC > StartTimeIndex_IC      %This ensures there is some data that exists to analyze
        
        Xlim1 = segmentTimeStart;
        Xlim2 = segmentTimeEnd;
        
         StartTimeIndex_IC = find(TimeAxis > Xlim1);
         StartTimeIndex_IC = StartTimeIndex_IC(1);

        EndTimeIndex_IC = find(TimeAxis > Xlim2);
        EndTimeIndex_IC = EndTimeIndex_IC(1);

            %%%%%Run the data through the Kalman Decoder Model
           CursorPositionOutput = PlayNeuralDataThroughKalmanEML(KM_features,StartTimeIndex_IC,EndTimeIndex_IC,A_KM,H_KM,K_KM);    

            Kalman_Data(:,StartTimeIndex_IC:EndTimeIndex_IC) = CursorPositionOutput;
  
    end  
   i  
end
    
save('ProcessedKalmanModelData.mat','Kalman_Data')
 
%% Try to look for replay:
%% Gather up all successful target trials and generate a peri-stimulus histogram of trajectories to get a mean trajectory template
%% Then look for that mean trajectory in the sleep data

%% Also generate templates for all of the "catch" trials and use those to generate a collection of "control" templates.

%% Cycle through each Simon Block creating a peri-stimulus spike histogram of successful trials based on the Kalman Model Output

SuccessfulTrialsKalmanData = struct;   %PeriStimulus data stored in this structure
SuccessCounter = 1;
SuccessfulCatchTrialsKalmanData = struct;
SuccessfulCatchTrialCounter = ones(20,1);

clear IndicesOfALLTrials
trialcounter = 1; 

for WhichSimonBlock = 1:10


    [GamePosition GameSequences GameClock GameTrialOutcomes GameNewTrialTimes GameNewTrialEndTimes Xlim1 Xlim2] = ExtractGameData(WhichSimonBlock);
    
 StartTimeIndex_IC = find(TimeAxis > Xlim1);
 StartTimeIndex_IC = StartTimeIndex_IC(1);

EndTimeIndex_IC = find(TimeAxis > Xlim2);
EndTimeIndex_IC = EndTimeIndex_IC(1)-1;

DateNum_IC_TimeAxis = datenum(TimeAxis);

%%%%%Run the data through the Kalman Decoder Model
CursorPositionOutput = PlayNeuralDataThroughKalmanEML(KM_features,StartTimeIndex_IC,EndTimeIndex_IC,A_KM,H_KM,K_KM);

figure(20010312)
clf
figure(20010312)
set(20010312,'Position',[50 100 1000 600])


subplot(4,1,1)

h = plotEEGSpectrogram(Xlim1,Xlim2,20010312,EEGSpectrogramTimeAxis,EEG_Spectrogram_Frequency_Axis,EEG_Spectrogram);

subplot(4,1,2)
S = plotRasterplot(StartTimeIndex_IC,EndTimeIndex_IC,20010312,DateNum_IC_TimeAxis,Thresholded_zSP);

subplot(4,1,3)

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
ylim([-0.4 0.4])
datetick('x','keeplimits');
ylabel('Cursor Position');
xlabel('Time');
title('Cursor Position')

subplot(4,1,4)

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
datetick('x','keeplimits');
ylabel('Magnitude');
xlabel('Time');
title('Raw Kalman Model Output')


Title_Details = sprintf(', %s to %s',datetime(Xlim1,'Format','HH-mm-ss'),datetime(Xlim2,'Format','HH-mm-ss'));
suptitle(['EEG, MEA, Cursor Position, and Kalman Filter Data' Title_Details])

PrintFigTitle = ['EEG, MEA, Cursor Position, and Kalman Filter Data' Title_Details  ' ' NeuralCodeSource '.pdf'];


fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 18 11];
fig.PaperSize = [18 11];
%print(PrintFigTitle,'-bestfit','-dpdf')


%%%Here is where the data to create the peri-stimulus activity plots is
%%%recorded:


for j = 1:length(GameNewTrialTimes)     %cycle through each trial on this block
    
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           IndicesOfALLTrials{trialcounter} = startIndex:endIndex;
            trialcounter = trialcounter+1;
   if GameTrialOutcomes(j) == 1         %if the trial was successful do 
       
       if sum((GameSequences(j,:) == [2 4 1 3])) == 4       %if the trial was the target sequence (which on this day was "2 4 1 3")
           
           %startIndex and endIndex are relative to the overall recording;
           %x_record(1) corresponds to the first time interval of the
           %normalized segment used for an individual round of the game so
           %the index (SIx and EIx) are normalized for each round of the
           %game
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
        SuccessfulTrialsKalmanData(SuccessCounter).KM_Data = CursorPositionOutput(:,SIx:EIx);
        SuccessfulTrialsKalmanData(SuccessCounter).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
        SuccessCounter = SuccessCounter+1;
       
       else
           
           
       end
       
       %%%These are the 20 different "Catch" trial sequences
       if sum((GameSequences(j,:) == [1 2 3 4])) == 4
           SequenceIndex = 1;
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).KM_Data = CursorPositionOutput(:,SIx:EIx);
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
            SuccessfulCatchTrialCounter(SequenceIndex) = SuccessfulCatchTrialCounter(SequenceIndex)+1;                
       end
       if sum((GameSequences(j,:) == [1 3 2 4])) == 4
           
           SequenceIndex = 2;
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).KM_Data = CursorPositionOutput(:,SIx:EIx);
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
            SuccessfulCatchTrialCounter(SequenceIndex) = SuccessfulCatchTrialCounter(SequenceIndex)+1;            
       end
       if sum((GameSequences(j,:) == [1 3 4 2])) == 4
           
           SequenceIndex = 3;
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).KM_Data = CursorPositionOutput(:,SIx:EIx);
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
            SuccessfulCatchTrialCounter(SequenceIndex) = SuccessfulCatchTrialCounter(SequenceIndex)+1;            
       end
       if sum((GameSequences(j,:) == [1 4 2 3])) == 4
           
           SequenceIndex = 4;
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).KM_Data = CursorPositionOutput(:,SIx:EIx);
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
            SuccessfulCatchTrialCounter(SequenceIndex) = SuccessfulCatchTrialCounter(SequenceIndex)+1;            
       end
       if sum((GameSequences(j,:) == [1 4 3 2])) == 4
           
           SequenceIndex = 5;
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).KM_Data = CursorPositionOutput(:,SIx:EIx);
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
            SuccessfulCatchTrialCounter(SequenceIndex) = SuccessfulCatchTrialCounter(SequenceIndex)+1;            
       end
       if sum((GameSequences(j,:) == [2 1 3 4])) == 4
           
           SequenceIndex = 6;
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).KM_Data = CursorPositionOutput(:,SIx:EIx);
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
            SuccessfulCatchTrialCounter(SequenceIndex) = SuccessfulCatchTrialCounter(SequenceIndex)+1;            
       end
       if sum((GameSequences(j,:) == [2 1 4 3])) == 4
           
           SequenceIndex = 7;
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).KM_Data = CursorPositionOutput(:,SIx:EIx);
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
            SuccessfulCatchTrialCounter(SequenceIndex) = SuccessfulCatchTrialCounter(SequenceIndex)+1;            
       end
       if sum((GameSequences(j,:) == [2 3 4 1])) == 4
           
           SequenceIndex = 8;
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).KM_Data = CursorPositionOutput(:,SIx:EIx);
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
            SuccessfulCatchTrialCounter(SequenceIndex) = SuccessfulCatchTrialCounter(SequenceIndex)+1;            
       end
       if sum((GameSequences(j,:) == [3 1 2 4])) == 4
           
           SequenceIndex = 9;
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).KM_Data = CursorPositionOutput(:,SIx:EIx);
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
            SuccessfulCatchTrialCounter(SequenceIndex) = SuccessfulCatchTrialCounter(SequenceIndex)+1;            
       end
       if sum((GameSequences(j,:) == [3 1 4 2])) == 4
           
           SequenceIndex = 10;
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).KM_Data = CursorPositionOutput(:,SIx:EIx);
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
            SuccessfulCatchTrialCounter(SequenceIndex) = SuccessfulCatchTrialCounter(SequenceIndex)+1;            
       end
       if sum((GameSequences(j,:) == [3 2 1 4])) == 4
           
           SequenceIndex = 11;
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).KM_Data = CursorPositionOutput(:,SIx:EIx);
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
            SuccessfulCatchTrialCounter(SequenceIndex) = SuccessfulCatchTrialCounter(SequenceIndex)+1;            
       end
        if sum((GameSequences(j,:) == [3 2 4 1])) == 4
           
           SequenceIndex = 12;
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).KM_Data = CursorPositionOutput(:,SIx:EIx);
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
            SuccessfulCatchTrialCounter(SequenceIndex) = SuccessfulCatchTrialCounter(SequenceIndex)+1;            
        end
        if sum((GameSequences(j,:) == [3 4 1 2])) == 4
           
           SequenceIndex = 13;
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).KM_Data = CursorPositionOutput(:,SIx:EIx);
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
            SuccessfulCatchTrialCounter(SequenceIndex) = SuccessfulCatchTrialCounter(SequenceIndex)+1;            
        end
        if sum((GameSequences(j,:) == [3 4 2 1])) == 4
           
           SequenceIndex = 14;
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).KM_Data = CursorPositionOutput(:,SIx:EIx);
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
            SuccessfulCatchTrialCounter(SequenceIndex) = SuccessfulCatchTrialCounter(SequenceIndex)+1;            
        end
        if sum((GameSequences(j,:) == [4 1 2 3])) == 4
           
           SequenceIndex = 15;
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).KM_Data = CursorPositionOutput(:,SIx:EIx);
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
            SuccessfulCatchTrialCounter(SequenceIndex) = SuccessfulCatchTrialCounter(SequenceIndex)+1;            
        end
        if sum((GameSequences(j,:) == [4 1 3 2])) == 4
           
           SequenceIndex = 16;
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).KM_Data = CursorPositionOutput(:,SIx:EIx);
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
            SuccessfulCatchTrialCounter(SequenceIndex) = SuccessfulCatchTrialCounter(SequenceIndex)+1;            
        end
        if sum((GameSequences(j,:) == [4 2 1 3])) == 4
           
           SequenceIndex = 17;
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).KM_Data = CursorPositionOutput(:,SIx:EIx);
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
            SuccessfulCatchTrialCounter(SequenceIndex) = SuccessfulCatchTrialCounter(SequenceIndex)+1;            
        end
        if sum((GameSequences(j,:) == [4 2 3 1])) == 4
           
           SequenceIndex = 18;
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).KM_Data = CursorPositionOutput(:,SIx:EIx);
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
            SuccessfulCatchTrialCounter(SequenceIndex) = SuccessfulCatchTrialCounter(SequenceIndex)+1;            
        end
        if sum((GameSequences(j,:) == [4 3 1 2])) == 4
           
           SequenceIndex = 19;
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).KM_Data = CursorPositionOutput(:,SIx:EIx);
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
            SuccessfulCatchTrialCounter(SequenceIndex) = SuccessfulCatchTrialCounter(SequenceIndex)+1;            
        end
        if sum((GameSequences(j,:) == [4 3 2 1])) == 4
           
           SequenceIndex = 20;
           [dummy startIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialTimes(j)));
           [dummy endIndex] = min(abs(DateNum_IC_TimeAxis - GameNewTrialEndTimes(j)));
          
           SIx = startIndex-StartTimeIndex_IC+1;
           EIx = endIndex-StartTimeIndex_IC+1;
           
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).KM_Data = CursorPositionOutput(:,SIx:EIx);
            SuccessfulCatchTrialsKalmanData(SequenceIndex).SuccessfulTrialsKD(SuccessfulCatchTrialCounter(SequenceIndex)).Timing = DateNum_IC_TimeAxis(startIndex:endIndex);
        
            SuccessfulCatchTrialCounter(SequenceIndex) = SuccessfulCatchTrialCounter(SequenceIndex)+1;            
       end
       
       
       
   end
  
end


end


%% Now assemble peri-stimulus histogram into single matrix and visualize the means of the KM dimensions, and use these to calculate the cross correlagrams with the simulated Kalman Model output

%%% First do the Target Sequence:
for go =1
% Determine which one is the longest and pad all other peri-stimulus firing
% curves with NaN to reach that length

for i = 1:length(SuccessfulTrialsKalmanData)
    
   spikeduration = SuccessfulTrialsKalmanData(i).Timing;
   duration_of_spikingData(i) = length(spikeduration); 
    
end
secondsStop = 4;
[longest_spike_train whichTrial] = max(duration_of_spikingData);

%Now gather all the data into a High D Matrix: N x M x P, where N = number
%of succseful trials, M = length of the activity vector, and P = number of
%channels (in this case dimensions of KM model):

AllTrialsKMData = zeros(length(SuccessfulTrialsKalmanData),3,longest_spike_train);

for i = 1:length(SuccessfulTrialsKalmanData)
    
   SpikeData = SuccessfulTrialsKalmanData(i).KM_Data;
    if length(SpikeData) < longest_spike_train
        SpikeData(:,end+1:longest_spike_train) = NaN; 
    end
    
   % size(SpikeData)
    AllTrialsKMData(i,:,:) = SpikeData;     %Dimensions are #trials X 3KM dimensions X duration
   
end

%%Find the mean and standard deviation of each KM dimension:

Mean_KM_One = nanmean(squeeze(AllTrialsKMData(:,1,:)));
Std_KM_One = nanstd(squeeze(AllTrialsKMData(:,1,:)));

Mean_KM_Two = nanmean(squeeze(AllTrialsKMData(:,2,:)));
Std_KM_Two = nanstd(squeeze(AllTrialsKMData(:,2,:)));

%Define the TimeAxis as a duration vector in seconds:

PSSH_TimeAxis = duration(seconds(0):seconds(0.02):seconds(0.02)*length(SuccessfulTrialsKalmanData(whichTrial).Timing));%,'format','hh:mm:ss');
[dummy endIndex] = min(abs(PSSH_TimeAxis - seconds(secondsStop)));


figure(100)
clf 
figure(100)
set(100,'Position',[50 350 1000 400])
set(100,'DefaultAxesFontSize',16);
%print(PrintFigTitle,'-dsvg','-painters')


subplot(1,2,1)
hold on
plot(PSSH_TimeAxis(1:endIndex),Mean_KM_One(1:endIndex),'b','LineWidth',2)
plot(PSSH_TimeAxis(1:endIndex),Mean_KM_One(1:endIndex)+Std_KM_One(1:endIndex)/sqrt(length(SuccessfulTrialsKalmanData)),'b','LineWidth',0.5)
plot(PSSH_TimeAxis(1:endIndex),Mean_KM_One(1:endIndex)-Std_KM_One(1:endIndex)/sqrt(length(SuccessfulTrialsKalmanData)),'b','LineWidth',0.5)
plot(PSSH_TimeAxis(1:endIndex),0*Mean_KM_One(1:endIndex),'r--')
xlabel('Time (seconds)')
xlim([seconds(0) seconds(secondsStop)])
title('KF X-Dimension Template')


subplot(1,2,2)
hold on
plot(PSSH_TimeAxis(1:endIndex),Mean_KM_Two(1:endIndex),'r','LineWidth',2)
plot(PSSH_TimeAxis(1:endIndex),Mean_KM_Two(1:endIndex)+Std_KM_Two(1:endIndex)/sqrt(length(SuccessfulTrialsKalmanData)),'r','LineWidth',0.5)
plot(PSSH_TimeAxis(1:endIndex),Mean_KM_Two(1:endIndex)-Std_KM_Two(1:endIndex)/sqrt(length(SuccessfulTrialsKalmanData)),'r','LineWidth',0.5)
plot(PSSH_TimeAxis(1:endIndex),0*Mean_KM_One(1:endIndex),'r--')
xlabel('Time (seconds)')
xlim([seconds(0) seconds(secondsStop)])
title('KF Y-Dimension Template')



PrintFigTitle = ['Mean + SE Kalman Dimensions During Successful Trial ' NeuralCodeSource '.svg'];
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 18 11];
% fig.PaperSize = [18 11];
% print(PrintFigTitle,'-bestfit','-dpdf')


print(PrintFigTitle,'-dsvg','-painters')

% Look for sliding coerllation coefficient over the the course of the recording including overnight:


[dummy endIndex] = min(abs(PSSH_TimeAxis - seconds(2)));

KM_One_Template = Mean_KM_One(1:endIndex);
KM_Two_Template = Mean_KM_Two(1:endIndex);

CrossCorr_KM1 = zeros(1,length(Kalman_Data(1,:)));
CrossCorr_KM2 = zeros(1,length(Kalman_Data(2,:)));

CrossCorr_KM1_pv = zeros(1,length(Kalman_Data(1,:)));
CrossCorr_KM2_pv = zeros(1,length(Kalman_Data(2,:)));


%Take the correlation of the 2 KMs and slide across the whole days
%and nights activity to find times that they are correlated:

%record both the raw correlation coefficient as well as the p-value of the
%correlation
% % 
tic
parfor i=1:(length(Kalman_Data(1,:)) - length(KM_One_Template))
    
    [r p] = corrcoef(KM_One_Template,Kalman_Data(1,i:i+length(KM_One_Template)-1));
    CrossCorr_KM1(i) = r(1,2);
    CrossCorr_KM1_pv(i) = p(1,2);
    
     [r p] = corrcoef(KM_Two_Template,Kalman_Data(2,i:i+length(KM_Two_Template)-1));
    CrossCorr_KM2(i) = r(1,2);
        CrossCorr_KM2_pv(i) = p(1,2);

    

end
toc

% save('CrossCorr_KMs_Templates.mat','CrossCorr_KM1','CrossCorr_KM2','CrossCorr_KM1_pv','CrossCorr_KM2_pv','KM_One_Template','KM_Two_Template')

end

%% Define the cross correlation thresholds for each KM1 dimension:

prctl_threshold=99;
CC_threshold_KM1 = prctile(CrossCorr_KM1,prctl_threshold);
CC_threshold_KM2 = prctile(CrossCorr_KM2,prctl_threshold);

%% Now find number of simultaneous threshold crossings during different Epochs

clear VeryInterestingIndices NumberVeryInterestIndices durationofcomparisonperiod totalfiringrateduringperiod

DateNum_IC_TimeAxis = datenum(TimeAxis);

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

for i = 1:6
    
    if i == 1
%The entire session
xLim1 = datetime(2020,11,23,23,25,0);
xLim2 = datetime(2020,11,25,9,0,0);
    end
    if i == 2
%%First night of sleep:
xLim1 = datetime(2020,11,23,23,30,0);
xLim2 = datetime(2020,11,24,9,0,0);
    end
    if i == 3
%The first pre-game rest period:
xLim1 = datetime(2020,11,24,13,50,0);
xLim2 = datetime(2020,11,24,14,25,0);
    end
    if i == 4
%Game playing:
xLim1 = datetime(2020,11,24,14,50,0);
xLim2 = datetime(2020,11,24,17,15,0);
    end
    if i == 5
%Post game play rest:
xLim1 = datetime(2020,11,24,17,15,0);
xLim2 = datetime(2020,11,24,17,50,0);
    end
    if i == 6
%Sleeping second overnight
xLim1 = datetime(2020,11,24,23,0,0);
xLim2 = datetime(2020,11,25,9,0,0);
    end
    

 TimeStartOfComparison = xLim1;
 TimeEndComparison = xLim2;

    StartTimeIndex_IC = find(TimeAxis > TimeStartOfComparison);
    StartTimeIndex_IC = StartTimeIndex_IC(1);

    EndTimeIndex_IC = find(TimeAxis > TimeEndComparison);
    EndTimeIndex_IC = EndTimeIndex_IC(1)-1;

    %%Periods of high correlation for forward target template

        InterestingIndices = find((CrossCorr_KM1(StartTimeIndex_IC:EndTimeIndex_IC) > CC_threshold_KM1) & (CrossCorr_KM2((StartTimeIndex_IC:EndTimeIndex_IC)) > CC_threshold_KM2));

        dInterestingIndices = diff(InterestingIndices);
        VeryInterestingIndices{i} = InterestingIndices(dInterestingIndices~=1);
        NumberVeryInterestIndices(i) = length(InterestingIndices(dInterestingIndices~=1));
        durationofcomparisonperiod(i) = hours(xLim2-xLim1);
        AllInterestingIndices = find((CrossCorr_KM1 > CC_threshold_KM1) & (CrossCorr_KM2 > CC_threshold_KM2));
        dInterestingIndices = diff(AllInterestingIndices);
        AllInterestingIndices = AllInterestingIndices(dInterestingIndices~=1);

  StartTimeIndex_IC = find(TimeAxis > xLim1);
StartTimeIndex_IC = StartTimeIndex_IC(1);

EndTimeIndex_IC = find(TimeAxis > xLim2);
EndTimeIndex_IC = EndTimeIndex_IC(1)-1;
  
totalfiringrateduringperiod(i) = sum(sum(SP(StartTimeIndex_IC:EndTimeIndex_IC,:)));
      
%% plot all of the instances of the cross correlation exceeding threshold
fignum = 2020+i;
figure(fignum)
close(fignum)
figure(fignum)

set(fignum,'Position',[50 -300 1200 1000])
set(fignum,'DefaultAxesFontSize',16);

StartTimeIndex_IC = find(TimeAxis > xLim1);
StartTimeIndex_IC = StartTimeIndex_IC(1);

EndTimeIndex_IC = find(TimeAxis > xLim2);
EndTimeIndex_IC = EndTimeIndex_IC(1)-1;


subplot(3,1,1)
h = plotEEGSpectrogram(xLim1,xLim2,fignum,EEGSpectrogramTimeAxis,EEG_Spectrogram_Frequency_Axis,EEG_Spectrogram);
t_lim = [datenum(xLim1) datenum(xLim2)];
xlim(t_lim);

pause(0.1)
% 
subplot(3,1,2)
%S = plotRasterplot(StartTimeIndex_IC,EndTimeIndex_IC,fignum,DateNum_IC_TimeAxis,Thresholded_zSP);
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),mean(SP(StartTimeIndex_IC:EndTimeIndex_IC,:),2))
t_lim = [datenum(xLim1) datenum(xLim2)];
xlim(t_lim);
datetick('x','keeplimits');
%ylabel('Mean Firing Rate')
pause(0.1)
% if i == 1 || i == 6
%     ylim([0 300])
% else
%     ylim([0 100])
% end

CorrelationThresholdCrossings = AllInterestingIndices;
CCrossings = DateNum_IC_TimeAxis*0;
CCrossings(CorrelationThresholdCrossings) = 2;
CCrossings = CCrossings+1;

subplot(3,1,3)
hold on
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),(CCrossings(StartTimeIndex_IC:EndTimeIndex_IC)))
%plot(EEGSpectrogramTimeAxis,(smootheddominantFrequency))
plot(EEGSpectrogramTimeAxis,(smoothedthresholdFrequency))
axis tight;
% ylim([0.1 12])
% if i == 6
%     ylim([0.1 7])
% end
t_lim = [datenum(xLim1) datenum(xLim2)];
%axis tight;
xlim(t_lim);

datetick('x','keeplimits');
ylabel('EEG 95% BandPower Freq');

Title_Details = sprintf(', %s to %s',datetime(xLim1,'Format','HH-mm-ss'),datetime(xLim2,'Format','HH-mm-ss'));
suptitle(['EEG, MEA, Dominant EEG Frequency, and CrossCorrelations' Title_Details])

PrintFigTitle = ['EEG, Mean Firing Rate, Dominant EEG Frequency, and CrossCorrelations' Title_Details  ' ' NeuralCodeSource 'PDF.pdf'];

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 24 24];
fig.PaperSize = [24 24];
%print(PrintFigTitle,'-dpdf')

%print(PrintFigTitle,'-dsvg','-painters')

end

%% Make some figures wth cross correlations:


%% A few during gameplay:
%% Cycle through each Block of the matching game, and record the the neural activity and timing of each trial
% taking special note of which of the trials were successfull vs. not and which were target vs. distractor


CorrelationThresholdCrossings = AllInterestingIndices;
CCrossings = DateNum_IC_TimeAxis*0;
CCrossings(CorrelationThresholdCrossings) = 3;
CCrossings = CCrossings-2;

for WhichSimonBlock = 1:10
    
    [GamePosition GameSequences GameClock GameTrialOutcomes GameNewTrialTimes GameNewTrialEndTimes Xlim1 Xlim2] = ExtractGameData(WhichSimonBlock)
    
 StartTimeIndex_IC = find(TimeAxis > Xlim1);
 StartTimeIndex_IC = StartTimeIndex_IC(1);

EndTimeIndex_IC = find(TimeAxis > Xlim2);
EndTimeIndex_IC = EndTimeIndex_IC(1);

DateNum_IC_TimeAxis = datenum(TimeAxis);

%%%%%Run the data through the Kalman Decoder Model
CursorPositionOutput = PlayNeuralDataThroughKalmanEML(KM_features,StartTimeIndex_IC,EndTimeIndex_IC,A_KM,H_KM,K_KM);


%%%Plot
figure(20010312)
clf
figure(20010312)
set(20010312,'Position',[50 -300 1600 2400])
set(20010312,'DefaultAxesFontSize',16);


t_lim = [datenum(Xlim1) datenum(Xlim2)];

 StartTimeIndex_IC = find(TimeAxis > Xlim1);
 StartTimeIndex_IC = StartTimeIndex_IC(1);

EndTimeIndex_IC = find(TimeAxis > Xlim2);
EndTimeIndex_IC = EndTimeIndex_IC(1)-1;

DateNum_IC_TimeAxis = datenum(TimeAxis);

%%%%%Run the data through the Kalman Decoder Model
CursorPositionOutput = PlayNeuralDataThroughKalmanEML(KM_features,StartTimeIndex_IC,EndTimeIndex_IC,A_KM,H_KM,K_KM);

subplot(5,1,1)

h = plotEEGSpectrogram(Xlim1,Xlim2,20010312,EEGSpectrogramTimeAxis,EEG_Spectrogram_Frequency_Axis,EEG_Spectrogram);

subplot(5,1,2)
S = plotRasterplot(StartTimeIndex_IC,EndTimeIndex_IC,20010312,DateNum_IC_TimeAxis,Thresholded_zSP);

subplot(5,1,3)

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
ylim([-0.4 0.4])
datetick('x','keeplimits');
ylabel('Cursor Position');
%xlabel('Time');
title('Cursor Position (Green: Successful Trial, Red: Unsuccessful Trial)')

subplot(5,1,4)

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
datetick('x','keeplimits');
ylabel('Position');
%xlabel('Time');
title('Kalman Filter Output (Light: Target Sequence, Dark: Distractor Sequence)')

% 

subplot(5,1,5)

hold on
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),CrossCorr_KM1(StartTimeIndex_IC:EndTimeIndex_IC),'b')
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),CrossCorr_KM2(StartTimeIndex_IC:EndTimeIndex_IC),'r')
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),0*CrossCorr_KM1(StartTimeIndex_IC:EndTimeIndex_IC) + CC_threshold_KM1,'b--')
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),0*CrossCorr_KM2(StartTimeIndex_IC:EndTimeIndex_IC) + CC_threshold_KM2,'r--')
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),CCrossings(StartTimeIndex_IC:EndTimeIndex_IC),'g--','linewidth',2)
axis tight;
xlim(t_lim)
ylim([-0.9 0.9])
datetick('x','keeplimits');
ylabel('Cross Correlation');
xlabel('Time');
title('Cross Correlation Between Templates and KF Output')

Title_Details = sprintf(', %s to %s',datetime(Xlim1,'Format','HH-mm-ss'),datetime(Xlim2,'Format','HH-mm-ss'));
%suptitle(['EEG, MEA, Cursor Position, and Kalman Filter Data' Title_Details])

PrintFigTitle = ['EEG, MEA, Cursor Position, Kalman Filter, and Correlation Data' Title_Details  '.svg'];
%PrintFigTitle = ['EEG, MEA, Cursor Position, Kalman Filter' Title_Details  '.pdf'];


% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 18 11];
% fig.PaperSize = [18 11];
% if printfigs == 1
% print(PrintFigTitle,'-bestfit','-dpdf')
% end

print(PrintFigTitle,'-dsvg','-painters')

end


%% And a few from overnight:

close all
CorrelationThresholdCrossings = AllInterestingIndices;
CCrossings = DateNum_IC_TimeAxis*0;
CCrossings(CorrelationThresholdCrossings) = 3;
CCrossings = CCrossings-2;

for WhichBlock = 1:4
    
    if WhichBlock == 1
        %random afternoony time
        Xlim1 = datetime(2020,11,24,13,55,0);
        Xlim2 = datetime(2020,11,24,13,58,0);
         
    end
    if WhichBlock == 2
        %random evening time
        Xlim1 = datetime(2020,11,24,18,10,0);
        Xlim2 = datetime(2020,11,24,18,13,0);  
    end
    if WhichBlock == 3
        %random overnight time1
        Xlim1 = datetime(2020,11,25,0,30,0);
        Xlim2 = datetime(2020,11,25,0,33,0);  
         
    end
    if WhichBlock == 4
        %random overnight time2
        Xlim1 = datetime(2020,11,25,0,43,0);
        Xlim2 = datetime(2020,11,25,0,46,0);  
        
         
    end
    
 StartTimeIndex_IC = find(TimeAxis > Xlim1);
 StartTimeIndex_IC = StartTimeIndex_IC(1);

EndTimeIndex_IC = find(TimeAxis > Xlim2);
EndTimeIndex_IC = EndTimeIndex_IC(1);

DateNum_IC_TimeAxis = datenum(TimeAxis);

%%%%%Run the data through the Kalman Decoder Model
CursorPositionOutput = PlayNeuralDataThroughKalmanEML(KM_features,StartTimeIndex_IC,EndTimeIndex_IC,A_KM,H_KM,K_KM);


%%%Plot
figure(200103)
clf
figure(200103)
set(200103,'Position',[50 0 1000 800])
t_lim = [datenum(Xlim1) datenum(Xlim2)];

subplot(4,1,1)
h = plotEEGSpectrogram(Xlim1,Xlim2,200103,EEGSpectrogramTimeAxis,EEG_Spectrogram_Frequency_Axis,EEG_Spectrogram);

subplot(4,1,2)
S = plotRasterplot(StartTimeIndex_IC,EndTimeIndex_IC,200103,DateNum_IC_TimeAxis,Thresholded_zSP);

subplot(4,1,3)

hold on
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),CursorPositionOutput(1:2,:))

axis tight;
xlim(t_lim)
%ylim([-220 220])
datetick('x','keeplimits');
ylabel('Magnitude');
xlabel('Time');
title('Raw Kalman Model Output')


subplot(4,1,4)

hold on
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),CrossCorr_KM1(StartTimeIndex_IC:EndTimeIndex_IC),'b')
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),CrossCorr_KM2(StartTimeIndex_IC:EndTimeIndex_IC),'r')
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),0*CrossCorr_KM1(StartTimeIndex_IC:EndTimeIndex_IC) + CC_threshold_KM1,'b--')
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),0*CrossCorr_KM2(StartTimeIndex_IC:EndTimeIndex_IC) + CC_threshold_KM2,'r--')
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),CCrossings(StartTimeIndex_IC:EndTimeIndex_IC),'g--','linewidth',2)
axis tight;
xlim(t_lim)
ylim([-0.9 0.9])
datetick('x','keeplimits');
ylabel('Cross Correlation');
xlabel('Time');


Title_Details = sprintf(', %s to %s',datetime(Xlim1,'Format','HH-mm-ss'),datetime(Xlim2,'Format','HH-mm-ss'));
suptitle(['EEG, MEA, Cursor Position, and Kalman Filter Data' Title_Details])

PrintFigTitle = ['EEG, MEA, Cursor Position, Kalman Filter, and Correlation Data' Title_Details  ' ' NeuralCodeSource '.pdf'];


fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 18 11];
fig.PaperSize = [18 11];
if printfigs == 1
print(PrintFigTitle,'-bestfit','-dpdf')
end

end

