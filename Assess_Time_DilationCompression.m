%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021

%% The goal of this script will be to open neural data, and assess for time dilation/compression by manipulating the parameters
% of the kalman filter model

%% Open intracranial neural data as well as surface EEG data:

clear all
close all
clc

printfigs = 0;       %Do (1) or do not (0) print figures as the analysis progresses

%% Add all directories to path to allow for ease of access to files
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

%Pad with zeros at discontinuities for plotting to prevent weird on the Raster plot:
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

load('ProcessedKalmanModelData.mat','Kalman_Data')
 
%% Workflow will be as follows: Generate multiple versions of the KD -- one at the original time constant (A), one at 0.75, one at 0.5 and one at 0.25 time constants
%% Then, for each time constant, generate a series of idealized templates: at logspace(log10(0.1),log10(10),21) relative speed using splines

%%
clear DistanceMetricsVector CorrelationMetricsVector

%Make a vector of compression factors:
%a = 1:10; b = fliplr(a); b = 1./b; b(end) = [];

b = [0.1 0.25 0.5 0.75 0.9];
a = [1.0 1.1 1.2 1.4 1.6 1.8 2 3 4 5 6 8 10];
templateCompressionFactor = [b a];

[DV, CVs] = CalculateDistanceFromTemplateOverTime(Kalman_Data,KM_One_Template,KM_Two_Template); %Serves as ground truth

KalmanData_adj = Kalman_Data; 

        for j = 1:length(templateCompressionFactor)
                        
            thisTemplateCompressionFactor = templateCompressionFactor(j);
            
            [KM_One_Template_adj KM_Two_Template_adj] = GenerateIdealizedTemplates(thisTemplateCompressionFactor,KM_One_Template,KM_Two_Template);  %Theoretical function to generate templates using splines
            
            [DistanceVector, CrossCorrVector] = CalculateDistanceFromTemplateOverTime(KalmanData_adj,KM_One_Template_adj,KM_Two_Template_adj);
            
            DistanceMetricsVector{j} = DistanceVector;
            CorrelationMetricsVector{j} = CrossCorrVector;
            
            j
            
        end


%save('DistanceVectors2D_ManyXChanges.mat','DistanceMetricsVector','CorrelationMetricsVector','DV','CVs','templateCompressionFactor','-v7.3')

%%
DistancePerecentileThreshold = 99;
  
CVS_threshold_KM1 = prctile(CVs(:,1),DistancePerecentileThreshold);
CVS_threshold_KM2 = prctile(CVs(:,2),DistancePerecentileThreshold);
 
Dis_threshold_KM1 = prctile(DV(:,1),DistancePerecentileThreshold);
Dis_threshold_KM2 = prctile(DV(:,2),DistancePerecentileThreshold);

%% Plot Some examples of distance vectors:
  
AllInterestingIndices = find((CVs(:,1) > CVS_threshold_KM1) & (CVs(:,2) > CVS_threshold_KM1));
dInterestingIndices = diff(AllInterestingIndices);
AllInterestingIndices = AllInterestingIndices(dInterestingIndices~=1);
        
CorrelationThresholdCrossings = AllInterestingIndices;
CCrossings = DateNum_IC_TimeAxis*0;
CCrossings(CorrelationThresholdCrossings) = 2;
CCrossings = CCrossings-1;

AllInterestingIndices = find((DV(:,1) > Dis_threshold_KM1) & (DV(:,2) > Dis_threshold_KM2));
dInterestingIndices = diff(AllInterestingIndices);
AllInterestingIndices = AllInterestingIndices(dInterestingIndices~=1);
        
DistanceThresholdCrossings = AllInterestingIndices;
DCrossings = DateNum_IC_TimeAxis*0;
DCrossings(DistanceThresholdCrossings) = 2;
DCrossings = DCrossings-1;        
  
for WhichSimonBlock = 2


    [GamePosition GameSequences GameClock GameTrialOutcomes GameNewTrialTimes GameNewTrialEndTimes Xlim1 Xlim2] = ExtractGameData(WhichSimonBlock);
    
 StartTimeIndex_IC = find(TimeAxis > Xlim1);
 StartTimeIndex_IC = StartTimeIndex_IC(1);

EndTimeIndex_IC = find(TimeAxis > Xlim2);
EndTimeIndex_IC = EndTimeIndex_IC(1)-1;

DateNum_IC_TimeAxis = datenum(TimeAxis);

%%%%%Run the data through the Kalman Decoder Model
CursorPositionOutput = PlayNeuralDataThroughKalmanEML(KM_features,StartTimeIndex_IC,EndTimeIndex_IC,A_KM,H_KM,K_KM);

fignum = 10+WhichSimonBlock;
figure(fignum)
clf
figure(fignum)
set(fignum,'Position',[50 100 1000 600])


subplot(4,1,1)

h = plotEEGSpectrogram(Xlim1,Xlim2,fignum,EEGSpectrogramTimeAxis,EEG_Spectrogram_Frequency_Axis,EEG_Spectrogram);

subplot(4,1,2)
S = plotRasterplot(StartTimeIndex_IC,EndTimeIndex_IC,fignum,DateNum_IC_TimeAxis,Thresholded_zSP);

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
   %pos(2) = -3;
   pos(3) = GameNewTrialEndTimes(j) - GameNewTrialTimes(j);
   pos(2) = -1;
   
   pos(4) = 2;  
   if sum((GameSequences(j,:) == [2 4 1 3])) == 4
       h = rectangle('Position',pos,'EdgeColor','none','FaceColor',[0.85 0.85 0.85]);
   else
       h = rectangle('Position',pos,'EdgeColor','none');
       h.FaceColor = [0.65 0.65 0.65];
   end
end

plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),DV(StartTimeIndex_IC:EndTimeIndex_IC,:))
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),CVs(StartTimeIndex_IC:EndTimeIndex_IC,:))
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),(CCrossings(StartTimeIndex_IC:EndTimeIndex_IC)),'LineWidth',10)
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),(DCrossings(StartTimeIndex_IC:EndTimeIndex_IC)),'LineWidth',2)
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
end



%%      
  clear VeryInterestingIndices NumberVeryInterestIndices durationofcomparisonperiod
  clear DistanceVeryInterestingIndices DistanceNumberVeryInterestIndices InterestingIndices_AdjustForTime
  clear DistanceAllInterestingIndicesModified CCAllInterestingIndicesModified
  
  DateNum_IC_TimeAxis = datenum(TimeAxis);
  
  tic
  
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

    %% Correlation Metric First

        InterestingIndices = find((CVs(StartTimeIndex_IC:EndTimeIndex_IC,1) > CVS_threshold_KM1) & (CVs(StartTimeIndex_IC:EndTimeIndex_IC,2) > CVS_threshold_KM2));  
        dInterestingIndices = diff(InterestingIndices);
        VeryInterestingIndices{i} = InterestingIndices(dInterestingIndices~=1);
        NumberVeryInterestIndices(i) = length(InterestingIndices(dInterestingIndices~=1));
        durationofcomparisonperiod(i) = hours(xLim2-xLim1);

        AllInterestingIndices = find((CVs(:,1) > CVS_threshold_KM1) & (CVs(:,2) > CVS_threshold_KM1));
        dInterestingIndices = diff(AllInterestingIndices);
        AllInterestingIndices = AllInterestingIndices(dInterestingIndices~=1);
        
        CorrelationThresholdCrossings = AllInterestingIndices;
        CCrossings = DateNum_IC_TimeAxis*0;
        CCrossings(CorrelationThresholdCrossings) = 2;
        CCrossings = CCrossings+1;

    %Now distance Metric:
        
        InterestingIndices = find((DV(StartTimeIndex_IC:EndTimeIndex_IC,1) > Dis_threshold_KM1) & (DV(StartTimeIndex_IC:EndTimeIndex_IC,2) > Dis_threshold_KM2));  
        dInterestingIndices = diff(InterestingIndices);
        DistanceVeryInterestingIndices{i} = InterestingIndices(dInterestingIndices~=1);
        DistanceNumberVeryInterestIndices(i) = length(InterestingIndices(dInterestingIndices~=1));
        
        AllInterestingIndices = find((DV(:,1) > Dis_threshold_KM1) & (DV(:,2) > Dis_threshold_KM2));
        dInterestingIndices = diff(AllInterestingIndices);
        DistanceAllInterestingIndices = AllInterestingIndices(dInterestingIndices~=1);
                    
        DistanceThresholdCrossings = DistanceAllInterestingIndices;
        DCrossings = DateNum_IC_TimeAxis*0;
        DCrossings(DistanceThresholdCrossings) = 2;
        DCrossings = DCrossings+1;    
    
    
        
  StartTimeIndex_IC = find(TimeAxis > xLim1);
StartTimeIndex_IC = StartTimeIndex_IC(1);

EndTimeIndex_IC = find(TimeAxis > xLim2);
EndTimeIndex_IC = EndTimeIndex_IC(1)-1;
  
totalfiringrateduringperiod(i) = sum(sum(SP(StartTimeIndex_IC:EndTimeIndex_IC,:)));
      
for j = 1:21
    
    modifyfactor = 1.0 + (j-1)*0.01;
    
    for k = 1:18
        
        %% Correlation Metric First

        CV_t = CorrelationMetricsVector{1,k};
        
        InterestingIndices = find((CV_t(StartTimeIndex_IC:EndTimeIndex_IC,1) > modifyfactor*CVS_threshold_KM1) & (CV_t(StartTimeIndex_IC:EndTimeIndex_IC,2) > modifyfactor*CVS_threshold_KM2));  
        dInterestingIndices = diff(InterestingIndices);
        
        InterestingIndices_AdjustForTime(j,k,i) = length(InterestingIndices(dInterestingIndices~=1));

        
                   
        AllInterestingIndices = find((CV_t(:,1) > modifyfactor*CVS_threshold_KM1) & (CV_t(:,2) > modifyfactor*CVS_threshold_KM2));
        dInterestingIndices = diff(AllInterestingIndices);
        CCAllInterestingIndicesModified{j,k} = AllInterestingIndices(dInterestingIndices~=1);
     
        %% Distance metric
        
        DV_t = DistanceMetricsVector{1,k};
        
        InterestingIndices = find((DV_t(StartTimeIndex_IC:EndTimeIndex_IC,1) > modifyfactor*Dis_threshold_KM1) & (DV_t(StartTimeIndex_IC:EndTimeIndex_IC,2) > modifyfactor*Dis_threshold_KM2));  
        dInterestingIndices = diff(InterestingIndices);
        
        DistanceInterestingIndices_AdjustForTime(j,k,i) = length(InterestingIndices(dInterestingIndices~=1));
                
        AllInterestingIndices = find((DV_t(:,1) > modifyfactor*Dis_threshold_KM1) & (DV_t(:,2) > modifyfactor*Dis_threshold_KM2));
        dInterestingIndices = diff(AllInterestingIndices);
        DistanceAllInterestingIndicesModified{j,k} = AllInterestingIndices(dInterestingIndices~=1);
        
        
    end
end

end

toc

%% Demonstrate some examples of activity at increased speed:

%Load the 40-125 Hz LFP band:

GammaBand = load('Band4_Array_2.mat');
GammaBand = GammaBand.Band4_array2;
DeltaBand = load('Band1_Array_2.mat');
DeltaBand = DeltaBand.Band1_array2;  


%%
%% Also analyze number of STCE at different quantifications of sleep stage:

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

STCE_EEGFreq_GamePlay = cell(18,1); STCE_EEGFreq_Sleep = cell(18,1);

for j = 1:18    %cycle through each speed
   
   CompressionFactor = templateCompressionFactor(j);
   TheseAreTheIndicesOfthe_STCEs = CCAllInterestingIndicesModified{1,j};   %these are indices of STCEs at this speed
   TimesOfAllSTCEs = TimeAxis(TheseAreTheIndicesOfthe_STCEs);
       
    %Game playing:
    xLim1 = datetime(2020,11,24,14,50,0);
    xLim2 = datetime(2020,11,24,17,15,0);
    
    a = find(TimesOfAllSTCEs > xLim1 & TimesOfAllSTCEs <xLim2);
    TimesOfGamePlayingSTCEs = TimesOfAllSTCEs(a);
    
    %Sleeping second overnight
    xLim1 = datetime(2020,11,24,23,0,0);
    xLim2 = datetime(2020,11,25,9,0,0);
    
    a = find(TimesOfAllSTCEs > xLim1 & TimesOfAllSTCEs <xLim2);
    TimesOfOvernightSTCEs = TimesOfAllSTCEs(a);
    
    EEGThresholdFrequencyAtThisSTCE = zeros(1,length(TimesOfGamePlayingSTCEs));
    for i = 1:length(TimesOfGamePlayingSTCEs)
        ThisSTCE = TimesOfGamePlayingSTCEs(i);
        [a TimeIndexOnEEGAxis] = min(abs(EEG_Spectrogram_DateTimeAxis - ThisSTCE));
        
        EEGThresholdFrequencyAtThisSTCE(i) = smoothedthresholdFrequency(TimeIndexOnEEGAxis);
    end
    STCE_EEGFreq_GamePlay{j} = EEGThresholdFrequencyAtThisSTCE;
    
    
    EEGThresholdFrequencyAtThisSTCE = zeros(1,length(TimesOfOvernightSTCEs));
    for i = 1:length(TimesOfOvernightSTCEs)
        ThisSTCE = TimesOfOvernightSTCEs(i);
        [a TimeIndexOnEEGAxis] = min(abs(EEG_Spectrogram_DateTimeAxis - ThisSTCE));
        
        EEGThresholdFrequencyAtThisSTCE(i) = smoothedthresholdFrequency(TimeIndexOnEEGAxis);
    end
    STCE_EEGFreq_Sleep{j} = EEGThresholdFrequencyAtThisSTCE;
    
end


% Look at overall frequency composition of each era as well:

xLim1 = datetime(2020,11,24,14,50,0);
xLim2 = datetime(2020,11,24,17,15,0); 
a = find(EEG_Spectrogram_DateTimeAxis > xLim1 & EEG_Spectrogram_DateTimeAxis <xLim2);
EEG_Frequencies_GameEpoch = smoothedthresholdFrequency(a);


xLim1 = datetime(2020,11,24,23,0,0);
xLim2 = datetime(2020,11,25,9,0,0);
a = find(EEG_Spectrogram_DateTimeAxis > xLim1 & EEG_Spectrogram_DateTimeAxis <xLim2);
EEG_Frequencies_Sleep = smoothedthresholdFrequency(a);

%% For my statistical test: will look at the total number of STCEs overnight. I will divide the overnight epoch into slow wave sleep and not slow wave sleep, and find the proportion of time in each. 
%My null hypothesis will be that the number of STCEs should be proportional
%to the time spent in each state (so if there are 85 STCEs and the
%overnight is 80% SWS and 20% not SWS, the null hypotehsis is 68 during SWS
%and 17 during NSWS

clear probsave SWSProp SWS_STCEsprop

speedIndices = [6 12 13 14];

for j = [1 2 3 4]
    
    STCEs = STCE_EEGFreq_Sleep{speedIndices(j)};
    
    TotalNumberSTCEs = length(STCEs);
    
   % clear SWSProp SWS_STCEsprop 
    cutoffFreq = 2:0.01:6;
    for i = 1:length(cutoffFreq)
        
        AllFrequencyCounts = length(EEG_Frequencies_Sleep);
        SWSFreqCounts = length(EEG_Frequencies_Sleep(EEG_Frequencies_Sleep<=cutoffFreq(i)));
        SWSP = SWSFreqCounts./AllFrequencyCounts;                
        SWSProp(j,i) = SWSFreqCounts./AllFrequencyCounts;
        SWS_STCEs = length(STCEs(STCEs<cutoffFreq(i)));
        SWS_STCEsprop(j,i) = length(STCEs(STCEs<cutoffFreq(i)))./TotalNumberSTCEs;
        
        y = 1 - binocdf(SWS_STCEs-1,TotalNumberSTCEs,SWSP);
        
        probsave(j,i) = y;
       
    end
    j
    TotalNumberSTCEs;
    SWS_STCEsprop;

% 
% figure(j)
% hold on
% plot(cutoffFreq,SWSProp)
% plot(cutoffFreq,SWS_STCEsprop)



%plot(cutoffFreq,probsave(j,:))
end

statisticalcutoff = find(probsave(1,:)<0.01);
statisticalcutoff = statisticalcutoff(1);
statisticalcutoffFrez = cutoffFreq(statisticalcutoff)


figure(29)
set(29,'DefaultAxesFontSize',16);
hold on
a = get(gca,'colororder');
plot(cutoffFreq,SWSProp(1,:),'k--')
plot(cutoffFreq,SWS_STCEsprop(1,:),'Color',a(1,:))
plot(cutoffFreq,SWS_STCEsprop(2,:),'Color',a(2,:))
plot(cutoffFreq,SWS_STCEsprop(3,:),'Color',a(3,:))
plot(cutoffFreq,SWS_STCEsprop(4,:),'Color',a(4,:))
legend('Slow Wave Sleep (SWS)','1x STCEs during SWS','2x STCEs during SWS','3x STCEs during SWS','4x STCEs during SWS','Location','SouthEast')
xlabel('SWS Cutoff Frequency')
ylabel('Proportion of Night')
ylim([0 1])



PrintFigTitle = ['Session 2 STCEs versus Proportion Of Night in SWS.svg'];
%  
% 
print(PrintFigTitle,'-dsvg','-painters')


%a = get(gca,'colororder');
%text(DateNum_IC_TimeAxis(StartTimeIndex_IC),3.5,'1x','FontSize',16,'Color',a(1,:))
%text(DateNum_IC_TimeAxis(StartTimeIndex_IC),5.5,'2x','FontSize',16,'Color',a(2,:))
%text(DateNum_IC_TimeAxis(StartTimeIndex_IC),7.5,'3x','FontSize',16,'Color',a(3,:))
%text(DateNum_IC_TimeAxis(StartTimeIndex_IC),9.5,'4x','FontSize',16,'Color',a(4,:))

%% plot all of the instances of the cross correlation exceeding threshold at a few speeds
%The entire session
xLim1 = datetime(2020,11,23,23,25,0);
xLim2 = datetime(2020,11,25,9,0,0);
%     end
%     if i == 2
% %%First night of sleep:
% xLim1 = datetime(2020,11,23,23,30,0);
% xLim2 = datetime(2020,11,24,9,0,0);
%     end
%     if i == 3
% %The first pre-game rest period:
% xLim1 = datetime(2020,11,24,13,50,0);
% xLim2 = datetime(2020,11,24,14,25,0);
%     end
%     if i == 4
% %Game playing:
% xLim1 = datetime(2020,11,24,14,50,0);
% xLim2 = datetime(2020,11,24,17,15,0);
%     end
%     if i == 5
% %Post game play rest:
% xLim1 = datetime(2020,11,24,17,15,0);
% xLim2 = datetime(2020,11,24,17,50,0);
%     end
%     if i == 6
%Sleeping second overnight
xLim1 = datetime(2020,11,24,23,0,0);
xLim2 = datetime(2020,11,25,9,0,0);
%     end
   

 TimeStartOfComparison = xLim1;
 TimeEndComparison = xLim2;

    StartTimeIndex_IC = find(TimeAxis > TimeStartOfComparison);
    StartTimeIndex_IC = StartTimeIndex_IC(1);

    EndTimeIndex_IC = find(TimeAxis > TimeEndComparison);
    EndTimeIndex_IC = EndTimeIndex_IC(1)-1;

    
speeds = [6 12 13 14]; 
CCrossingsX = zeros(4,length(DateNum_IC_TimeAxis));

for z = 1:4
    TheseAreTheIndicesOfthe_STCEs = CCAllInterestingIndicesModified{1,speeds(z)};
    
    CCrossings = DateNum_IC_TimeAxis*0;
    CCrossings(TheseAreTheIndicesOfthe_STCEs) = 1;
    CCrossings = CCrossings+z*2;
    
    CCrossingsX(z,:) = CCrossings;
end

%Define which times are above or below a threshold for SWS:

cutoffFreq = 3;
% 


SWSMarker = smoothedthresholdFrequency < cutoffFreq;
SWSMarker = double(SWSMarker);

fignum = 2020;
figure(fignum)
set(fignum,'Position',[50 100 1200 600])
set(fignum,'DefaultAxesFontSize',16);

subplot(2,1,1)
h = plotEEGSpectrogram(xLim1,xLim2,fignum,EEGSpectrogramTimeAxis,EEG_Spectrogram_Frequency_Axis,EEG_Spectrogram);
t_lim = [datenum(xLim1) datenum(xLim2)];
xlim(t_lim);

subplot(2,1,2)



hold on

plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),CCrossingsX(:,StartTimeIndex_IC:EndTimeIndex_IC))
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),0*CCrossingsX(1,StartTimeIndex_IC:EndTimeIndex_IC)+2,'Color',[1 1 1],'LineWidth',2)
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),0*CCrossingsX(1,StartTimeIndex_IC:EndTimeIndex_IC)+4,'Color',[1 1 1],'LineWidth',2)
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),0*CCrossingsX(1,StartTimeIndex_IC:EndTimeIndex_IC)+6,'Color',[1 1 1],'LineWidth',2)
plot(DateNum_IC_TimeAxis(StartTimeIndex_IC:EndTimeIndex_IC),0*CCrossingsX(1,StartTimeIndex_IC:EndTimeIndex_IC)+8,'Color',[1 1 1],'LineWidth',2)
plot(EEGSpectrogramTimeAxis,(smoothedthresholdFrequency))
%plot(EEGSpectrogramTimeAxis,SWSMarker+1,'k')
a = get(gca,'colororder');
text(DateNum_IC_TimeAxis(StartTimeIndex_IC),3.5,'1x','FontSize',16,'Color',a(1,:))
text(DateNum_IC_TimeAxis(StartTimeIndex_IC),5.5,'2x','FontSize',16,'Color',a(2,:))
text(DateNum_IC_TimeAxis(StartTimeIndex_IC),7.5,'3x','FontSize',16,'Color',a(3,:))
text(DateNum_IC_TimeAxis(StartTimeIndex_IC),9.5,'4x','FontSize',16,'Color',a(4,:))
axis tight;



t_lim = [datenum(xLim1) datenum(xLim2)];
xlim(t_lim);
ylim([0 10])

datetick('x','keeplimits');
ylabel('EEG 95% BandPower Freq');
 

PrintFigTitle = ['STCEs versus 95 percent band power at 4 speeds.svg'];
%  
% 
%print(PrintFigTitle,'-dsvg','-painters')
% 
