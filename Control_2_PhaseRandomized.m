%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021

%here we take the output of the kalman decoder and generate new random
%kalman filter outputs with the same frequency content (so it has the same
%kind of peaks and valleys) but different phase (so it breaks any non-random
%temporal structure like replay of specific sequences).

%%

load('AnalysisWorkspace.mat')

kalmanOutput = Kalman_Data(1:2,:)';
kalmanOutput(isnan(kalmanOutput)) = 0;

%everything is easier if there are an even number of data points
if mod(size(kalmanOutput,1),2)==1
    kalmanOutput = kalmanOutput(1:(end-1),:);
end

%%
nRepeats = 100;

clear VeryInterestingIndices_randomphase NumberVeryInterestIndices_randomphase durationofcomparisonperiod_randomphase

for z = 1:nRepeats

%randomize phase by changing the complex angle of the FFT output but
%preserving the magnitude
fftOut = fft(kalmanOutput);
randTheta = rand((size(fftOut,1)/2)-1,2)*2*pi;
randComplex = cos(randTheta) + sin(randTheta)*1i;

newFFT = zeros(size(fftOut));
newFFT(2:(end/2),:) = abs(fftOut(2:(end/2),:)).*randComplex;
newFFT(((end/2)+2):end,:) = abs(fftOut(((end/2)+2):end,:)).*flipud(conj(randComplex));

randTimeSeries = ifft(newFFT);

% %%
% %compare random time series to original
% figure
% subplot(2,1,1);
% hold on;
% plot(kalmanOutput(:,1));
% plot(randTimeSeries(:,1));
% xlabel('Time Step');
% ylabel('X Kalman Output');
% legend({'Original','Synthetic'});
% 
% subplot(2,1,2);
% hold on;
% plot(kalmanOutput(:,2));
% plot(randTimeSeries(:,2));
% xlabel('Time Step');
% ylabel('Y Kalman Output');

%%

randTimeSeries;

Kalman_Data_RandomPhase = randTimeSeries';

CrossCorr_KM1_randomphase = zeros(1,length(Kalman_Data_RandomPhase(1,:)));
CrossCorr_KM2_randomphase = zeros(1,length(Kalman_Data_RandomPhase(2,:)));

CrossCorr_KM1_pv_randomphase = zeros(1,length(Kalman_Data_RandomPhase(1,:)));
CrossCorr_KM2_pv_randomphase = zeros(1,length(Kalman_Data_RandomPhase(2,:)));

%Take the correlation of the KMs and slide across the whole days
%and nights activity to find times that they are correlated:

%record both the raw correlation coefficient as well as the p-value of the
%correlation
% 
tic
parfor i=1:(length(Kalman_Data_RandomPhase(1,:)) - length(KM_One_Template))
    
    [r p] = corrcoef(KM_One_Template,Kalman_Data_RandomPhase(1,i:i+length(KM_One_Template)-1));
    CrossCorr_KM1_randomphase(i) = r(1,2);
    CrossCorr_KM1_pv_randomphase(i) = p(1,2);
    
     [r p] = corrcoef(KM_Two_Template,Kalman_Data_RandomPhase(2,i:i+length(KM_Two_Template)-1));
    CrossCorr_KM2_randomphase(i) = r(1,2);
    CrossCorr_KM2_pv_randomphase(i) = p(1,2);

end
toc

%Use the same thresholds as forward:
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

        InterestingIndices = find((CrossCorr_KM1_randomphase(StartTimeIndex_IC:EndTimeIndex_IC) > CC_threshold_KM1_random) & (CrossCorr_KM2_randomphase((StartTimeIndex_IC:EndTimeIndex_IC)) > CC_threshold_KM2_random));

        dInterestingIndices = diff(InterestingIndices);
        VeryInterestingIndices_randomphase{z,i} = InterestingIndices(dInterestingIndices~=1);
        NumberVeryInterestIndices_randomphase(z,i) = length(InterestingIndices(dInterestingIndices~=1));
        durationofcomparisonperiod_randomphase(z,i) = hours(xLim2-xLim1);

        AllInterestingIndices = find((CrossCorr_KM1_randomphase > CC_threshold_KM1_random) & (CrossCorr_KM2_randomphase > CC_threshold_KM2_random));
        dInterestingIndices = diff(AllInterestingIndices);
        AllInterestingIndices = AllInterestingIndices(dInterestingIndices~=1);
end



z

save('OverlapsPhaseRandomization_100repeats.mat','NumberVeryInterestIndices_randomphase')

end

code1 = now;
save(['OverlapsPhaseRandomization_100repeats' num2str(code1) '.mat'],'NumberVeryInterestIndices_randomphase')

%%

%% Stats comparing random to real data:
RealDataNumberAt99Percent_night1 = 4;
RealDataNumberAt99Percent_night2 = 85;

NumberTCEsGeneratedByRandomPhaseControl_night1 = NumberVeryInterestIndices_randomphase(:,2);
NumberTCEsGeneratedByRandomPhaseControl_night2 = NumberVeryInterestIndices_randomphase(:,6);

%%% Compare Both Nights at the same time:
[pranksum_night1 h stats1] = ranksum(NumberTCEsGeneratedByRandomPhaseControl_night1,RealDataNumberAt99Percent_night1,'method','exact');
[pranksum_night2 h stats2] = ranksum(NumberTCEsGeneratedByRandomPhaseControl_night2,RealDataNumberAt99Percent_night2,'method','exact');
[pranksum_night1vs2 h stats12] = ranksum(NumberTCEsGeneratedByRandomPhaseControl_night1,NumberTCEsGeneratedByRandomPhaseControl_night2);%,'method','exact');

[pranksum_night1 h stats1] = ranksum(NumberTCEsGeneratedByRandomPhaseControl_night1,RealDataNumberAt99Percent_night1,'method','exact');
[pranksum_night2 h stats2] = ranksum(NumberTCEsGeneratedByRandomPhaseControl_night2,RealDataNumberAt99Percent_night2,'method','exact');
[pranksum_night1vs2 h stats12] = ranksum(NumberTCEsGeneratedByRandomPhaseControl_night1,NumberTCEsGeneratedByRandomPhaseControl_night2);%,'method','exact');


figure(2123)

set(2123,'Position',[150 350 1200 450])


clf
hold on
set(gca,'FontSize',12)
h1 = histogram(NumberTCEsGeneratedByRandomPhaseControl_night1,'BinMethod','integers')
h2 = histogram(NumberTCEsGeneratedByRandomPhaseControl_night2,'BinMethod','integers')
line([RealDataNumberAt99Percent_night1 RealDataNumberAt99Percent_night1],[0 max(max(h1.Values),max(h2.Values))],'LineWidth',1,'Color','b','LineStyle','-.')
line([RealDataNumberAt99Percent_night2 RealDataNumberAt99Percent_night2],[0 max(max(h1.Values),max(h2.Values))],'LineWidth',1,'Color','r','LineStyle','-.')
text(20,10,['Night 1 Distribution vs 4 TCEs, Ranksum test, p = ' num2str(pranksum_night1)],'FontSize',12)
text(20,8,['Night 2 Distribution vs 85 TCEs, Ranksum test, p = ' num2str(pranksum_night2) '*'],'FontSize',12)
text(20,6,['Night 1 Distribution vs Night 2 Distribution, p = ' num2str(pranksum_night1vs2)],'FontSize',12)
text(20,14','Night 1 (Before): 4 Target TCEs','Color','b','FontSize',12)
text(60,14','Night 2 (After): 85 Target TCEs','Color','r','FontSize',12)

ylim([0 max(max(h1.Values),max(h2.Values))])

xlabel('Number of Threshold Crossing Events (TCEs) Observed')
ylabel('Frequency (Bin Count)')
legend('Night 1 Bootstrap Bin Counts','Night 2 Bootstrap Bin Counts','Location','SouthEast')

PrintFigTitle = ['Threshold Crossing Events Using Phase-Randomized Shuffled Data Night 2' num2str(code1) '.pdf'];

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 14 12];
fig.PaperSize = [14 12];

print(PrintFigTitle,'-bestfit','-dpdf')

%%
%% Repeat the phase randomization process above, but limiting the randomization to X minute long segments of activity
%%
%%
nRepeats = 100;
randomization_time_window = 30;     %how many minute chunks to use for phase randomization scheme

clear VeryInterestingIndices_randomphase NumberVeryInterestIndices_randomphase durationofcomparisonperiod_randomphase

for z = 1:nRepeats

tic
    clear STI ETI
randTimeSeries = [];

SessionStart = datetime(2020,11,23,23,0,0);
SessionEnd = datetime(2020,11,25,9,0,0);
DurationTotalRecording = SessionEnd-SessionStart;
Number_of_TimeSegments = ceil(minutes(DurationTotalRecording)/randomization_time_window);

    for y = 1:Number_of_TimeSegments
    SegmentStart = SessionStart + minutes((y-1)*randomization_time_window);
    SegmentEnd = SegmentStart + minutes(randomization_time_window);
    if SegmentEnd > SessionEnd
       SegmentEnd = SessionEnd;
    end

    StartTimeIndex_IC = find(TimeAxis > SegmentStart);
    StartTimeIndex_IC = StartTimeIndex_IC(1);

%     STI(y) = StartTimeIndex_IC;
    
    EndTimeIndex_IC = find(TimeAxis > SegmentEnd);
    EndTimeIndex_IC = EndTimeIndex_IC(1)-1;
    
%     ETI(y) = EndTimeIndex_IC
%     y
%     EndTimeIndex_IC - StartTimeIndex_IC
%     
    if EndTimeIndex_IC > StartTimeIndex_IC
        
            
        %randomize phase by changing the complex angle of the FFT output but
        %preserving the magnitude
        kalmanOutput_this_segment = kalmanOutput(StartTimeIndex_IC:EndTimeIndex_IC,:);

        %everything is easier if there are an even number of data points
        if mod(size(kalmanOutput_this_segment,1),2)==1
            kalmanOutput_this_segment = kalmanOutput_this_segment(1:(end-1),:);
        end

        fftOut = fft(kalmanOutput_this_segment);
        randTheta = rand((size(fftOut,1)/2)-1,2)*2*pi;
        randComplex = cos(randTheta) + sin(randTheta)*1i;

        newFFT = zeros(size(fftOut));
        newFFT(2:(end/2),:) = abs(fftOut(2:(end/2),:)).*randComplex;
        newFFT(((end/2)+2):end,:) = abs(fftOut(((end/2)+2):end,:)).*flipud(conj(randComplex));

        randTimeSeries_this_segment = ifft(newFFT);

        if isempty(randTimeSeries)
            randTimeSeries = randTimeSeries_this_segment;
        else
            randTimeSeries = [randTimeSeries;randTimeSeries_this_segment];
        end
        
    end
    
    end

%%
toc
'made the data. now do the correlation'

Kalman_Data_RandomPhase = randTimeSeries';

CrossCorr_KM1_randomphase = zeros(1,length(Kalman_Data_RandomPhase(1,:)));
CrossCorr_KM2_randomphase = zeros(1,length(Kalman_Data_RandomPhase(2,:)));

CrossCorr_KM1_pv_randomphase = zeros(1,length(Kalman_Data_RandomPhase(1,:)));
CrossCorr_KM2_pv_randomphase = zeros(1,length(Kalman_Data_RandomPhase(2,:)));

%Take the correlation of the KMs and slide across the whole days
%and nights activity to find times that they are correlated:

%record both the raw correlation coefficient as well as the p-value of the
%correlation
% 

parfor i=1:(length(Kalman_Data_RandomPhase(1,:)) - length(KM_One_Template))
    
    [r p] = corrcoef(KM_One_Template,Kalman_Data_RandomPhase(1,i:i+length(KM_One_Template)-1));
    CrossCorr_KM1_randomphase(i) = r(1,2);
    CrossCorr_KM1_pv_randomphase(i) = p(1,2);
    
     [r p] = corrcoef(KM_Two_Template,Kalman_Data_RandomPhase(2,i:i+length(KM_Two_Template)-1));
    CrossCorr_KM2_randomphase(i) = r(1,2);
    CrossCorr_KM2_pv_randomphase(i) = p(1,2);

end
toc

%Use the same thresholds as forward:
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

    if EndTimeIndex_IC > length(CrossCorr_KM1_randomphase)
        EndTimeIndex_IC = length(CrossCorr_KM1_randomphase);
    end
    
    %%Periods of high correlation for forward target template

        InterestingIndices = find((CrossCorr_KM1_randomphase(StartTimeIndex_IC:EndTimeIndex_IC) > CC_threshold_KM1_random) & (CrossCorr_KM2_randomphase((StartTimeIndex_IC:EndTimeIndex_IC)) > CC_threshold_KM2_random));

        dInterestingIndices = diff(InterestingIndices);
        VeryInterestingIndices_randomphase{z,i} = InterestingIndices(dInterestingIndices~=1);
        NumberVeryInterestIndices_randomphase(z,i) = length(InterestingIndices(dInterestingIndices~=1));
        durationofcomparisonperiod_randomphase(z,i) = hours(xLim2-xLim1);

        AllInterestingIndices = find((CrossCorr_KM1_randomphase > CC_threshold_KM1_random) & (CrossCorr_KM2_randomphase > CC_threshold_KM2_random));
        dInterestingIndices = diff(AllInterestingIndices);
        AllInterestingIndices = AllInterestingIndices(dInterestingIndices~=1);
end

z

save(['OverlapsPhaseRandomization_' num2str(randomization_time_window) 'minutes_window_.mat'],'NumberVeryInterestIndices_randomphase')

end

code1 = now;
save(['OverlapsPhaseRandomization_100repeats_' num2str(randomization_time_window) 'minutes_window_' num2str(code1) '.mat'],'NumberVeryInterestIndices_randomphase')
%%
%% Stats comparing random to real data:
RealDataNumberAt99Percent_night1 = 4;
RealDataNumberAt99Percent_night2 = 85;

NumberTCEsGeneratedByRandomPhaseControl_night1 = NumberVeryInterestIndices_randomphase(:,2);
NumberTCEsGeneratedByRandomPhaseControl_night2 = NumberVeryInterestIndices_randomphase(:,6);



%%% Compare Both Nights at the same time:
[pranksum_night1 h stats1] = ranksum(NumberTCEsGeneratedByRandomPhaseControl_night1,RealDataNumberAt99Percent_night1,'method','exact');
[pranksum_night2 h stats2] = ranksum(NumberTCEsGeneratedByRandomPhaseControl_night2,RealDataNumberAt99Percent_night2,'method','exact');
[pranksum_night1vs2 h stats12] = ranksum(NumberTCEsGeneratedByRandomPhaseControl_night1,NumberTCEsGeneratedByRandomPhaseControl_night2);%,'method','exact');

[h pttest_night1] = ttest2(NumberTCEsGeneratedByRandomPhaseControl_night1,RealDataNumberAt99Percent_night1);
[h pttest_night2] = ttest2(NumberTCEsGeneratedByRandomPhaseControl_night2,RealDataNumberAt99Percent_night2);
[h pttest_night1vs2] = ttest2(NumberTCEsGeneratedByRandomPhaseControl_night1,NumberTCEsGeneratedByRandomPhaseControl_night2);%,'method','exact');

n1 = NumberTCEsGeneratedByRandomPhaseControl_night1;
n2 = NumberTCEsGeneratedByRandomPhaseControl_night2;

figure(2123)

set(2123,'Position',[150 350 1200 450])
set(2123,'DefaultAxesFontSize',16);

clf
hold on
set(gca,'FontSize',16)
h1 = histogram(NumberTCEsGeneratedByRandomPhaseControl_night1,'BinMethod','integers')
h2 = histogram(NumberTCEsGeneratedByRandomPhaseControl_night2,'BinMethod','integers')
line([RealDataNumberAt99Percent_night1 RealDataNumberAt99Percent_night1],[0 max(max(h1.Values),max(h2.Values))],'LineWidth',1,'Color','b','LineStyle','-.')
line([RealDataNumberAt99Percent_night2 RealDataNumberAt99Percent_night2],[0 max(max(h1.Values),max(h2.Values))],'LineWidth',1,'Color','r','LineStyle','-.')

text(25,6,['Night 1 Distribution vs 4 TCEs, Ranksum, p = ' num2str(pranksum_night1)],'FontSize',16)
text(25,5,['Night 2 Distribution vs 85 TCEs, Ranksum, p = ' num2str(pranksum_night2) '*'],'FontSize',16)
text(25,4,['Night 1 Distribution vs Night 2 Distribution, p = ' num2str(pranksum_night1vs2) '*'],'FontSize',16)

text(25,12,['Night 1 Median Random TCEs: ' num2str(median(n1)) ' (' num2str(prctile(n1,25)) ' - '  num2str(prctile(n1,75)) ')'],'FontSize',16)
text(25,11,['Night 2 Median Random TCEs: ' num2str(median(n2)) ' (' num2str(prctile(n2,25)) ' - '  num2str(prctile(n2,75)) ')'],'FontSize',16)

text(25,10,['Night 1 Distribution vs 4 TCEs, t test, p = ' num2str(pttest_night1)],'FontSize',16)
text(25,9,['Night 2 Distribution vs 85 TCEs, t test, p = ' num2str(pttest_night2) '*'],'FontSize',16)
text(25,8,['Night 1 Distribution vs Night 2 Distribution, p = ' num2str(pttest_night1vs2) '*'],'FontSize',16)


text(15,14','Night 1 (Before): 4 Target TCEs','Color','b','FontSize',16)
text(50,14','Night 2 (After): 85 Target TCEs','Color','r','FontSize',16)

ylim([0 max(max(h1.Values),max(h2.Values))])

xlabel('Number of Threshold Crossing Events (TCEs) Observed')
ylabel('Frequency (Bin Count)')
legend('Night 1 Bootstrap Bin Counts','Night 2 Bootstrap Bin Counts','Location','SouthEast')

PrintFigTitle = ['Threshold Crossing Events Using Phase-Randomized Shuffled Data Night 2_' num2str(now) 'minutes_window.svg'];

%fig = gcf;
%fig.PaperUnits = 'inches';
%fig.PaperPosition = [0 0 14 12];
%fig.PaperSize = [14 12];

%print(PrintFigTitle,'-bestfit','-dpdf')

print(PrintFigTitle,'-dsvg','-painters')
%set(2123,'DefaultAxesFontSize',16);