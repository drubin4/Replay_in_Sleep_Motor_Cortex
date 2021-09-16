%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021

%% How often do we expect a random motif to appear?

%Compare to controls where segments were sped up or slowed down by various
%amounts

load('AnalysisWorkspace.mat')

KM_One_Template_Target = KM_One_Template;
KM_Two_Template_Target = KM_Two_Template;

b = [0.1 0.25 0.5 0.75 0.9];
a = [1.0 1.1 1.2 1.4 1.6 1.8 2 3 4 5 6 8 10];
templateCompressionFactor = [b a];


for hijk = 1:length(templateCompressionFactor)

CompressionFactor = templateCompressionFactor(hijk);
% %% Take 100 different random segments from the pre-game rest of the same duration as the template
% %% and see how often they appear overnight to generate a control for statistical comparison

%create EEG 95% frequency band value for plotting:
for do = 1
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
end


clear VeryInterestingIndices_random NumberVeryInterestIndices_random durationofcomparisonperiod_random
nRepeats = 100;
AllRandomCrossCor1 = cell(nRepeats,1);
AllRandomCrossCor2 = cell(nRepeats,1);

for z = 1:nRepeats

% Identify a random segment from the 30 minute pre-game nap time
% (which has a mix of asleep and awake EEG activity)

tNapStart = datetime(2020,11,24,13,53,0);
tNapEnd = datetime(2020,11,24,14,22,0);

StartTimeIndex_IC = find(TimeAxis > tNapStart);
StartTimeIndex_IC = StartTimeIndex_IC(1);

EndTimeIndex_IC = find(TimeAxis > tNapEnd);
EndTimeIndex_IC = EndTimeIndex_IC(1)-1;

KalmanDataDuringNap = Kalman_Data(:,StartTimeIndex_IC:EndTimeIndex_IC);

duration_to_randomly_select = length(KM_One_Template);

Random_start_point = randi(length(KalmanDataDuringNap)-duration_to_randomly_select-1);

KM_One_Template_Random = KalmanDataDuringNap(1,Random_start_point:Random_start_point+duration_to_randomly_select-1);
KM_Two_Template_Random = KalmanDataDuringNap(2,Random_start_point:Random_start_point+duration_to_randomly_select-1);

[KM_One_Template_adj, KM_Two_Template_adj] = GenerateIdealizedTemplates(CompressionFactor,KM_One_Template_Random,KM_Two_Template_Random);  %Theoretical function to generate templates using splines
           
CrossCorr_KM1_random = zeros(1,length(Kalman_Data(1,:)));
CrossCorr_KM2_random = zeros(1,length(Kalman_Data(2,:)));


%record both the raw correlation coefficient as well as the p-value of the
%correlation
% 
tic
parfor i=1:(length(Kalman_Data(1,:)) - length(KM_One_Template_adj))
    
    [r p] = corrcoef(KM_One_Template_adj,Kalman_Data(1,i:i+length(KM_One_Template_adj)-1));
    CrossCorr_KM1_random(i) = r(1,2);
    
     [r p] = corrcoef(KM_Two_Template_adj,Kalman_Data(2,i:i+length(KM_Two_Template_adj)-1));
    CrossCorr_KM2_random(i) = r(1,2);   
    
end
toc

AllRandomCrossCor1{z} = CrossCorr_KM1_random;
AllRandomCrossCor2{z} = CrossCorr_KM2_random;



CC_threshold_KM1_random = CC_threshold_KM1;
CC_threshold_KM2_random = CC_threshold_KM2;

%% Now find number of simultaneous threshold crossings during different Epochs

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

        InterestingIndices = find((CrossCorr_KM1_random(StartTimeIndex_IC:EndTimeIndex_IC) > CC_threshold_KM1_random) & (CrossCorr_KM2_random((StartTimeIndex_IC:EndTimeIndex_IC)) > CC_threshold_KM2_random));

        dInterestingIndices = diff(InterestingIndices);
        VeryInterestingIndices_random{z,i} = InterestingIndices(dInterestingIndices~=1);
        NumberVeryInterestIndices_random(z,i) = length(InterestingIndices(dInterestingIndices~=1));
        durationofcomparisonperiod_random(z,i) = hours(xLim2-xLim1);

        AllInterestingIndices = find((CrossCorr_KM1_random > CC_threshold_KM1_random) & (CrossCorr_KM2_random > CC_threshold_KM2_random));
        dInterestingIndices = diff(AllInterestingIndices);
        AllInterestingIndices = AllInterestingIndices(dInterestingIndices~=1);
end

z

end

fname = ['WithCounts_OverlapsRandomSegment_100repeats_CompressionFactor_' num2str(CompressionFactor) 'x_' num2str(now) '.mat']
save(fname,'NumberVeryInterestIndices_random','NumberVeryInterestIndices_random','AllRandomCrossCor1','AllRandomCrossCor2','-v7.3');

%%
%% Stats comparing random to real data:
load('TimeDilationCompressionSTCECounts.mat');

code1 = now;

RealDataNumberAt99Percent_night1 = InterestingIndices_AdjustForTime(1,hijk,2);
RealDataNumberAt99Percent_night2 = InterestingIndices_AdjustForTime(1,hijk,6);

rn1 = RealDataNumberAt99Percent_night1;
rn2 = RealDataNumberAt99Percent_night2;

NumberTCEsGeneratedByRandomSegmentsat99Percent_night1 = NumberVeryInterestIndices_random(:,2);
NumberTCEsGeneratedByRandomSegmentsat99Percent_night2 = NumberVeryInterestIndices_random(:,6);


%%% Compare Both Nights at the same time:
[pranksum_night1 h stats1] = ranksum(NumberTCEsGeneratedByRandomSegmentsat99Percent_night1,RealDataNumberAt99Percent_night1,'method','exact');
[pranksum_night2 h stats2] = ranksum(NumberTCEsGeneratedByRandomSegmentsat99Percent_night2,RealDataNumberAt99Percent_night2,'method','exact');
[pranksum_night1vs2 h stats12] = ranksum(NumberTCEsGeneratedByRandomSegmentsat99Percent_night1,NumberTCEsGeneratedByRandomSegmentsat99Percent_night2,'method','exact');

[h, pttest_night1] = ttest2(NumberTCEsGeneratedByRandomSegmentsat99Percent_night1,RealDataNumberAt99Percent_night1);
[h, pttest_night2] = ttest2(NumberTCEsGeneratedByRandomSegmentsat99Percent_night2,RealDataNumberAt99Percent_night2);
[h, pttest_night1vs2] = ttest2(NumberTCEsGeneratedByRandomSegmentsat99Percent_night1,NumberTCEsGeneratedByRandomSegmentsat99Percent_night2);


figure(2123)

set(2123,'Position',[150 350 1200 450])
set(2123,'DefaultAxesFontSize',16);

n1 = NumberTCEsGeneratedByRandomSegmentsat99Percent_night1;
n2 = NumberTCEsGeneratedByRandomSegmentsat99Percent_night2;

clf
hold on
set(gca,'FontSize',16)
h1 = histogram(NumberTCEsGeneratedByRandomSegmentsat99Percent_night1,'BinMethod','integers')
h2 = histogram(NumberTCEsGeneratedByRandomSegmentsat99Percent_night2,'BinMethod','integers')
line([RealDataNumberAt99Percent_night1 RealDataNumberAt99Percent_night1],[0 max(max(h1.Values),max(h2.Values))],'LineWidth',1,'Color','b','LineStyle','-.')
line([RealDataNumberAt99Percent_night2 RealDataNumberAt99Percent_night2],[0 max(max(h1.Values),max(h2.Values))],'LineWidth',1,'Color','r','LineStyle','-.')


text(5,40,['Night 1 Median Random TCEs: ' num2str(median(n1)) ' (' num2str(prctile(n1,25)) ' - '  num2str(prctile(n1,75)) ')'],'FontSize',16)
text(5,35,['Night 2 Median Random TCEs: ' num2str(median(n2)) ' (' num2str(prctile(n2,25)) ' - '  num2str(prctile(n2,75)) ')'],'FontSize',16)
text(5,30,['Night 1 Distribution vs ' num2str(rn1) ' TCEs, Ranksum test, p = ' num2str(pranksum_night1)],'FontSize',16)
text(5,25,['Night 2 Distribution vs ' num2str(rn2) ' TCEs, Ranksum test, p = ' num2str(pranksum_night2) '*'],'FontSize',16)
text(5,20,['Night 1 Distribution vs Night 2 Distribution, p = ' num2str(pranksum_night1vs2)],'FontSize',16)

text(5,70,['Night 1 Distribution vs ' num2str(rn1) ' TCEs, T test, p = ' num2str(pttest_night1)],'FontSize',16)
text(5,65,['Night 2 Distribution vs ' num2str(rn2) ' TCEs, T test, p = ' num2str(pttest_night2) '*'],'FontSize',16)
text(5,60,['Night 1 Distribution vs Night 2 Distribution, p = ' num2str(pttest_night1vs2)],'FontSize',16)


text(5,85',['Night 1 (Before): ' num2str(rn1) ' Target TCEs'],'Color','b','FontSize',16)
text(50,85',['Night 2 (After): ' num2str(rn2) ' Target TCEs'],'Color','r','FontSize',16)

ylim([0 max(max(h1.Values),max(h2.Values))])
%xlim([-1 (max([n1; n2; rn1; rn2])+1)])

xlabel('Number of Threshold Crossing Events (TCEs) Observed')
ylabel('Frequency (Bin Count)')
legend('Night 1 Bootstrap Bin Counts','Night 2 Bootstrap Bin Counts','Location','Best')





PrintFigTitle = ['Number of TCEs from Randomly Sampled Epochs_' num2str(CompressionFactor) 'x_' num2str(now) '.svg'];

print(PrintFigTitle,'-dsvg','-painters')


end

%% Look at all the generated data
S = what;
S=  S.mat;

counter = 1;
clear Night1 Night2 DataName

for i = 1:length(S)
   
    fname = S(i);
    tic
    if length(fname{1}) > 62
        if strcmp('WithCounts_OverlapsRandomSegment_100repeats_CompressionFactor_',fname{1}(1:62))
           
            load(fname{1});
            NumberTCEsGeneratedByRandomSegmentsat99Percent_night1 = NumberVeryInterestIndices_random(:,2);
            NumberTCEsGeneratedByRandomSegmentsat99Percent_night2 = NumberVeryInterestIndices_random(:,6);

            Night1{counter} = NumberTCEsGeneratedByRandomSegmentsat99Percent_night1;
            Night2{counter} = NumberTCEsGeneratedByRandomSegmentsat99Percent_night2;
            DataName{counter} = fname{1};
            
            counter = counter+1
        end
    end
    
    toc
end

%%
clear mean_n1 median_n1 fiveprct_n1 twentyfiveprct_n1 seventyfiveprct_n1 ninetyfiveprct_n1
clear mean_n2 median_n2 fiveprct_n2 twentyfiveprct_n2 seventyfiveprct_n2 ninetyfiveprct_n2
clear CompressionSize

for i = 1:length(Night1)
    FileName = DataName{i};
   FileName(1:62) = [];
   x = find(FileName == 'x');
   FileName(x:end) = [];
   
   CompressionSize(i) = str2num(FileName);  
   
end

[SortexCSize ind] = unique(CompressionSize);
ind(12) = 17
ind(13) = 19

for i = 1:length(ind)
   mean_n1(i) = mean(Night1{ind(i)}); 
   std_n1(i) = std(Night1{ind(i)});
   median_n1(i) = median(Night1{ind(i)});
   fiveprct_n1(i) = prctile(Night1{ind(i)},5);
   twentyfiveprct_n1(i) = prctile(Night1{ind(i)},25);
   seventyfiveprct_n1(i) = prctile(Night1{ind(i)},75);
   ninetyfiveprct_n1(i) = prctile(Night1{ind(i)},95);
   min_n1(i) = min(Night1{ind(i)});
   max_n1(i) = max(Night1{ind(i)});
   
   mean_n2(i) = mean(Night2{ind(i)}); 
   std_n2(i) = std(Night2{ind(i)});
   median_n2(i) = median(Night2{ind(i)});
   fiveprct_n2(i) = prctile(Night2{ind(i)},5);
   twentyfiveprct_n2(i) = prctile(Night2{ind(i)},25);
   seventyfiveprct_n2(i) = prctile(Night2{ind(i)},75);
   ninetyfiveprct_n2(i) = prctile(Night2{ind(i)},95);
         
   min_n2(i) = min(Night2{ind(i)});
   max_n2(i) = max(Night2{ind(i)});
end

load('TimeDilationCompressionSTCECounts.mat');

RealDataNumberAt99Percent_night1 = InterestingIndices_AdjustForTime(1,:,2);
RealDataNumberAt99Percent_night2 = InterestingIndices_AdjustForTime(1,:,6);

for i = 1:length(RealDataNumberAt99Percent_night2)
       
    [h pttest] = ttest2(Night2{ind(i)},RealDataNumberAt99Percent_night2(i));
    [pranksun h] = ranksum(Night2{ind(i)},RealDataNumberAt99Percent_night2(i),'method','exact');

    TTest_p_n2(i) = pttest;
    RankSum_p_n2(i) = pranksun
    
end


fignum = 27;
figure(fignum)
set(fignum,'Position',[150 100 1200 600])
set(fignum,'DefaultAxesFontSize',16);

ninetyfivepercentCI = mean_n2 + 1.96*std_n2./sqrt(length(Night2{1}));
fivepercentCI = mean_n2 - 1.96*std_n2./sqrt(length(Night2{1}));

hold on
plot(SortexCSize,RealDataNumberAt99Percent_night2,'-o')
plot(SortexCSize,median_n2,'-*')
plot(SortexCSize,seventyfiveprct_n2,'-d')
plot(SortexCSize,ninetyfiveprct_n2,'-^')
plot(SortexCSize,1./RankSum_p_n2,'-p')
plot(SortexCSize,(SortexCSize*0)+(1/0.05),'--')
%plot(SortexCSize,mean_n1,'-+')
%plot(SortexCSize,fiveprct_n1,'-x')
%plot(SortexCSize,twentyfiveprct_n1,'-s')

%plot(SortexCSize,ninetyfivepercentCI,'-p')

legend('Real Data','Median (control)','75th% (control)','95th% (control)','1/p-value (RankSum)','Location','SouthEast')
title('Replay vs. Time Compression, Compared to Random Control Distribution')
ylabel('Number of Replay Events')
xlabel('Time Compression Factor')
%set(gca,'YScale','log')


PrintFigTitle = ['Night 2 Replay vs. Time Compression, Compared to Random Control Distribution Log.svg'];

%print(PrintFigTitle,'-dsvg','-painters')
%%

clear mean_n1 median_n1 fiveprct_n1 twentyfiveprct_n1 seventyfiveprct_n1 ninetyfiveprct_n1
clear mean_n2 median_n2 fiveprct_n2 twentyfiveprct_n2 seventyfiveprct_n2 ninetyfiveprct_n2
clear CompressionSize

for i = 1:length(Night1)
    FileName = DataName{i};
   FileName(1:62) = [];
   x = find(FileName == 'x');
   FileName(x:end) = [];
   
   CompressionSize(i) = str2num(FileName);  
   
end

[SortexCSize ind] = unique(CompressionSize);
ind(12) = 17
ind(13) = 19

for i = 1:length(ind)
   mean_n1(i) = mean(Night1{ind(i)}); 
   std_n1(i) = std(Night1{ind(i)});
   median_n1(i) = median(Night1{ind(i)});
   fiveprct_n1(i) = prctile(Night1{ind(i)},5);
   twentyfiveprct_n1(i) = prctile(Night1{ind(i)},25);
   seventyfiveprct_n1(i) = prctile(Night1{ind(i)},75);
   ninetyfiveprct_n1(i) = prctile(Night1{ind(i)},95);
   
   mean_n2(i) = mean(Night2{ind(i)}); 
   std_n2(i) = std(Night2{ind(i)});
   median_n2(i) = median(Night2{ind(i)});
   fiveprct_n2(i) = prctile(Night2{ind(i)},5);
   twentyfiveprct_n2(i) = prctile(Night2{ind(i)},25);
   seventyfiveprct_n2(i) = prctile(Night2{ind(i)},75);
   ninetyfiveprct_n2(i) = prctile(Night2{ind(i)},95);
      
end

load('TimeDilationCompressionSTCECounts.mat');

RealDataNumberAt99Percent_night1 = InterestingIndices_AdjustForTime(1,:,2);
RealDataNumberAt99Percent_night2 = InterestingIndices_AdjustForTime(1,:,6);

for i = 1:length(RealDataNumberAt99Percent_night2)
       
    [h pttest] = ttest2(Night1{ind(i)},RealDataNumberAt99Percent_night1(i));
    [pranksun h] = ranksum(Night1{ind(i)},RealDataNumberAt99Percent_night1(i),'method','exact');

    TTest_p_n1(i) = pttest;
    RankSum_p_n1(i) = pranksun
    
end


fignum = 27;
figure(fignum)
set(fignum,'Position',[150 100 1200 600])
set(fignum,'DefaultAxesFontSize',16);

ninetyfivepercentCI = mean_n1 + 1.96*std_n1./sqrt(length(Night1{1}));
fivepercentCI = mean_n1 - 1.96*std_n1./sqrt(length(Night1{1}));

hold on
plot(SortexCSize,RealDataNumberAt99Percent_night1,'-o')
plot(SortexCSize,median_n1,'-*')
plot(SortexCSize,seventyfiveprct_n1,'-d')
plot(SortexCSize,ninetyfiveprct_n1,'-^')
plot(SortexCSize,1./RankSum_p_n1,'-p')
plot(SortexCSize,(SortexCSize*0)+(1/0.05),'--')
%plot(SortexCSize,mean_n1,'-+')
%plot(SortexCSize,fiveprct_n1,'-x')
%plot(SortexCSize,twentyfiveprct_n1,'-s')

%plot(SortexCSize,ninetyfivepercentCI,'-p')

legend('Real Data','Median (control)','75th% (control)','95th% (control)','1/p-value (RankSum)','Location','SouthEast')
title('Replay vs. Time Compression, Compared to Random Control Distribution')
ylabel('Number of Replay Events')
xlabel('Time Compression Factor')
%set(gca,'YScale','log')


PrintFigTitle = ['Night 1 Replay vs. Time Compression, Compared to Random Control Distribution.svg'];

%print(PrintFigTitle,'-dsvg','-painters')
%%

fignum = 28;
figure(fignum)
set(fignum,'Position',[150 100 1200 600])
set(fignum,'DefaultAxesFontSize',16);

hold on
plot(SortexCSize,RealDataNumberAt99Percent_night2,'-o','LineWidth',2)
%plot(SortexCSize,min_n2,'k--','MarkerFaceColor','k')
plot(SortexCSize,fiveprct_n2,'k--','MarkerFaceColor','k')
plot(SortexCSize,twentyfiveprct_n2,'b-','MarkerFaceColor','b')
plot(SortexCSize,median_n2,'r-','MarkerFaceColor','r')
plot(SortexCSize,seventyfiveprct_n2,'b-','MarkerFaceColor','b')
%plot(SortexCSize,max_n2,'k--','MarkerFaceColor','k')
plot(SortexCSize,ninetyfiveprct_n2,'k--','MarkerFaceColor','k')

for i = 1:length(SortexCSize)
    SortexCSize_X = randn(100,1)*0.05 +SortexCSize(i);
    plot(SortexCSize_X,Night2{ind(i)},'.')
    
end

legend('Real Data','5%/05% (control)','25%/75% (control)','Median (control)','Location','NorthEast')
title('Replay vs. Time Compression, Compared to Phase Randomized Control Distribution')
ylabel('Number of Replay Events')
xlabel('Time Compression Factor')
set(gca,'YScale','log')
xlim([0 10.1])
ylim([1 round(1.1*max(RealDataNumberAt99Percent_night2))])


PrintFigTitle = ['Night 2 Replay vs. Time Compression, Compared to All  Random Control Data Distribution Log.svg'];

print(PrintFigTitle,'-dsvg','-painters')


fignum = 29;
figure(fignum)
set(fignum,'Position',[150 100 1200 600])
set(fignum,'DefaultAxesFontSize',16);

hold on
plot(SortexCSize,RealDataNumberAt99Percent_night1,'-o','LineWidth',2)
%plot(SortexCSize,min_n1,'k--','MarkerFaceColor','k')
plot(SortexCSize,fiveprct_n1,'k--','MarkerFaceColor','k')
plot(SortexCSize,twentyfiveprct_n1,'b-','MarkerFaceColor','b')
plot(SortexCSize,median_n1,'r-','MarkerFaceColor','r')
plot(SortexCSize,seventyfiveprct_n1,'b-','MarkerFaceColor','b')
%plot(SortexCSize,max_n1,'k--','MarkerFaceColor','k')
plot(SortexCSize,ninetyfiveprct_n2,'k--','MarkerFaceColor','k')

for i = 1:length(SortexCSize)
    SortexCSize_X = randn(100,1)*0.05 +SortexCSize(i);
    plot(SortexCSize_X,Night1{ind(i)},'.')
    
end

legend('Real Data','5%/95% (control)','25%/75% (control)','Median (control)','Location','NorthEast')
title('Replay vs. Time Compression, Compared to Phase Randomized Control Distribution')
ylabel('Number of Replay Events')
xlabel('Time Compression Factor')
set(gca,'YScale','log')
xlim([0 10.1])
ylim([1 round(1.1*max(RealDataNumberAt99Percent_night1))])

PrintFigTitle = ['Night 1 Replay vs. Time Compression, Compared to All Random Control Data Distribution Log.svg'];

print(PrintFigTitle,'-dsvg','-painters')