%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021

%% The Goal of This Script will be to Open each and every NS5 file, and attempt to find the time of each ripple event
%% I will use the criteria outlined in Jiang et al. (2020) hippocampus to start
%% First part will be to find the time of each Ripple. If succesfful, may consider also looking for the time of each Sharp Wave Riple
clear all
close all
clc

%% Add all directories to path to allow for ease of access to NS5 files

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
    
    
%% Generate a list of ALL the NS5 files and store it in cell array "NS5_filenames": first for lateral array, then medial
clear NS5_filenames_lateral

NS5dir = '*DataDirectory*\Data\_Lateral\NSP Data';
NSFiles = dir(NS5dir);

j = 1;
for i = 1:length(NSFiles)
    putativefilename = NSFiles(i).name;
    
    
    if length(putativefilename) > 2
        if strcmp(putativefilename(end-2:end),'ns5')
            NS5_filenames_lateral{j} = putativefilename;
            j = j+1;
        end
    end
end

clear NS5_filenames_medial

NS5dir = '*DataDirectory*\Data\_Medial\NSP Data';
NSFiles = dir(NS5dir);

j = 1;
for i = 1:length(NSFiles)
    putativefilename = NSFiles(i).name;
    
    
    if length(putativefilename) > 2
        if strcmp(putativefilename(end-2:end),'ns5')
            NS5_filenames_medial{j} = putativefilename;
            j = j+1;
        end
    end
end

%% Open each NS5 file

for exampleScriptnumber = 1:length(NS5_filenames_medial)

openthisNS5fileM = NS5_filenames_medial{exampleScriptnumber};
openthisNS5fileL = NS5_filenames_lateral{exampleScriptnumber};

NS5dirM = '*DataDirectory*\Data\_Medial\NSP Data\';
NS5dirL = '*DataDirectory*\Data\_Lateral\NSP Data\';

filepathM = [NS5dirM openthisNS5fileM];
filepathL = [NS5dirL openthisNS5fileL];

sM = dir(filepathM);
sL = dir(filepathL);

if (sM.bytes < 4e9) && (sL.bytes < 4e9) && (sM.bytes > 1e3) && (sL.bytes > 1e3)

%ReadNSx_Seconds(filepath)
A = openNSx(filepathM,'uv');
B = openNSx(filepathL,'uv');
    
%% Steps to find ripples:

%1) Downsample to 1000 Hz
%2) Bandpass 60-120 Hz (sixth-order Butterworth IIR band-pass filter, with zero-phaseforward and reverse filtering)
%3) Root mean square (RMS) over time of the filtered signal using a moving average of 20 ms
%4) 80th percentile of RMS values for each channel ia set as heuristic cut-off -- whenever a channel's signal exceeds this cut-off, a putative ripple event is detected
%5) Adjacent putative ripple event indices lessthan 40 ms apart were merged, with the center of the new event chosenby the highest RMS value within the merged time window
%6) Each putative ripple event was then evaluated based on the number of distinct peaks in the LFP signal (low-passed at 120 Hz) surrounding the event center
%-- A 40-ms time bin was shifted (5 ms per shift) across ±50 ms, and at least one such time bin must include more than three peaks 
%-- The first and the last peak cannot be both less than 7 ms away from the edges of this time bin for the putative ripple event to be considered for subsequentanalyses. 
%-- In addition, the distance between two consecutive ripple centers must exceed 40 ms. 
%7) To determine the duration of each ripple, eft edge was defined as the first time point before the ripple centerwhere RMS value exceeded the percentile threshold used for initial rip-ple detection, and the right edge was defined as the first time point afterthe ripple center where RMS value fell below the threshold; all othertime points within the duration of this ripple must exceed the RMSthreshold.
%%

MedialArrayData = A.Data{2};
LateralArrayData = B.Data{2};

MedialArrayData(97,:) = [];
LateralArrayData(97,:) = [];


%Define start times from the NS5 data file
MedialStartTime = A.MetaTags.DateTimeRaw;
MedialStartTime(5) = MedialStartTime(5)-5;
MedialStartTime(3) = [];
MedialStartTime(6) = MedialStartTime(6)+MedialStartTime(7)/1000;
MedialStartTime(7) = [];
MStartTime = datetime(MedialStartTime);

LateralStartTime = B.MetaTags.DateTimeRaw;
LateralStartTime(5) = LateralStartTime(5)-5;
LateralStartTime(3) = [];
LateralStartTime(6) = LateralStartTime(6)+LateralStartTime(7)/1000;
LateralStartTime(7) = [];
LStartTime = datetime(LateralStartTime);

MedialTimeSync = size(A.Data{1},2)./A.MetaTags.SamplingFreq;
LateralTimeSync = size(B.Data{1},2)./B.MetaTags.SamplingFreq;

UnifiedStartTimeM = MStartTime + seconds(MedialTimeSync);
UnifiedStartTimeL = LStartTime + seconds(LateralTimeSync);

lengthofMedialTimeAxis = size(MedialArrayData,2)-1;
TimeAxisM = UnifiedStartTimeM + seconds((0:1:lengthofMedialTimeAxis)*(1/A.MetaTags.SamplingFreq));

lengthofLateralTimeAxis = size(LateralArrayData,2)-1;
TimeAxisL = UnifiedStartTimeL + seconds((0:1:lengthofLateralTimeAxis)*(1/B.MetaTags.SamplingFreq));

%% 1) Downsample everything to 1000 Hz:

MedialArrayData_ds = downsample(MedialArrayData',(A.MetaTags.SamplingFreq)/1000)';
LateralArrayData_ds = downsample(LateralArrayData',(B.MetaTags.SamplingFreq)/1000)';
TimeAxisM_ds = downsample(TimeAxisM,(A.MetaTags.SamplingFreq)/1000);
TimeAxisL_ds = downsample(TimeAxisL,(B.MetaTags.SamplingFreq)/1000);

%% %2) Bandpass 60-120 Hz (sixth-order Butterworth IIR band-pass filter, with zero-phaseforward and reverse filtering)

Filtered_MedialArrayData_ds = zeros(size(MedialArrayData_ds));
Filtered_LateralArrayData_ds = zeros(size(LateralArrayData_ds));

LP_Filtered_MedialArrayData_ds = zeros(size(MedialArrayData_ds));
LP_Filtered_LateralArrayData_ds = zeros(size(LateralArrayData_ds));

tic
for i = 1:96
    Filtered_MedialArrayData_ds(i,:) = bandpass(MedialArrayData_ds(i,:),[60 100],1000,'ImpulseResponse','iir');
    Filtered_LateralArrayData_ds(i,:) = bandpass(LateralArrayData_ds(i,:),[60 100],1000,'ImpulseResponse','iir');
    
    LP_Filtered_MedialArrayData_ds(i,:) = lowpass(MedialArrayData_ds(i,:),120,1000,'ImpulseResponse','iir');
    LP_Filtered_LateralArrayData_ds(i,:) = lowpass(LateralArrayData_ds(i,:),120,1000,'ImpulseResponse','iir');
    
end
toc

%% 3) Root mean square (RMS) over time of the filtered signal using a moving average of 20 ms
Filtered_MedialArrayData_rms = zeros(size(MedialArrayData_ds));
Filtered_LateralArrayData_rms = zeros(size(LateralArrayData_ds));

for i = 1:96
    A = Filtered_MedialArrayData_ds(i,:);
    Filtered_MedialArrayData_rms(i,:) = sqrt(movsum(A.^2,[0 20])/20);
   
    A = Filtered_LateralArrayData_ds(i,:);
    Filtered_LateralArrayData_rms(i,:) = sqrt(movsum(A.^2,[0 20])/20);
end

%% %4) 80th percentile of RMS values for each channel ia set as heuristic cut-off -- whenever a channel's signal exceeds this cut-off, a putative ripple event is detected
%5) Adjacent putative ripple event indices lessthan 40 ms apart were merged, with the center of the new event chosenby the highest RMS value within the merged time window
Filtered_MedialArrayData_RippleYesNo = zeros(size(MedialArrayData_ds));
MedialArrayRipples = struct;

for i = 1:96

    cutoff = prctile(Filtered_MedialArrayData_rms(i,:),90);     %anywhere the filtered RMS is >90% is putative Ripple
    
    pRipple = find(Filtered_MedialArrayData_rms(i,:)>cutoff);   %These are all the indices of "putative ripple"
    dRipple = diff(pRipple);                                    %this is the number of indices (= milliseconds because 1000 Hz sampling) between Ripple indices
    dRipple = [dRipple 1];
       
    RippleIndices = cell([sum(dRipple>40)+1 1]);
    rippleID = 1;

    
    for j = 1:length(pRipple)
        
        if dRipple(j) > 40                                                      %if the next ripple is more tha 40ms apart, it's a new ripple
        rippleID = rippleID+1;
        end
        RippleIndices{rippleID}  = [RippleIndices{rippleID}, pRipple(j)];       %Add indices onto this Ripple
          
    end
    
    %6) Each putative ripple event was then evaluated based on the number of distinct peaks in the LFP signal (low-passed at 120 Hz) surrounding the event center
    %-- A 40-ms time bin was shifted (5 ms per shift) across ±50 ms, and at least one such time bin must include more than three peaks 
    %-- The first and the last peak cannot be both less than 7 ms away from the edges of this time bin for the putative ripple event to be considered for subsequentanalyses. 
    %-- In addition, the distance between two consecutive ripple centers must exceed 40 ms. 
    %7) To determine the duration of each ripple, eft edge was defined as the first time point before the ripple centerwhere RMS value exceeded the percentile threshold used for initial rip-ple detection, and the right edge was defined as the first time point afterthe ripple center where RMS value fell below the threshold; all othertime points within the duration of this ripple must exceed the RMSthreshold.
    realRippleCounter = 0;
    RealRipples = [];
    
    for k = 1:length(RippleIndices)
        
        putativeRippleindices = RippleIndices{k};
        
        putatitveRippleLP120Hz = LP_Filtered_MedialArrayData_ds(i,putativeRippleindices);
        
        putativeRippleisARealRipple = 0;
        %Need to check for 3 peaks within 40 ms; to have three peaks        
        %a peak will be defined as where the diff of the sign of the diff
        %of the waveform is -2
        %Additionally, if either the first and last peak need to be more than
        %7 ms from the edge of the time bin, then the ripple must be at least 10 ms long 

        %Only consider segments that are at least 41 ms long:
        
        if length(putatitveRippleLP120Hz) > 40        
            for m = 1:length(putatitveRippleLP120Hz)-41
               dPR = diff(putatitveRippleLP120Hz(m:m+40));
               signchanges = diff(sign(dPR));
               peaks = find(signchanges == -2);     %signchanges == -2 are the indices where a peak occurs (where the first differential goes from positive to negative)
               if length(peaks) >2
                   if peaks(1) >7 || peaks(3) < (length(signchanges) - 7)
                        putativeRippleisARealRipple = 1;
                   end
               end
            end
        end    
          
        %%%if I wanted to also consider shorter ripples I could include this
        %code:
%         else
%            if length(putatitveRippleLP120Hz) > 10      
%                dPR = diff(putatitveRippleLP120Hz);            
%                signchanges = diff(sign(dPR));
%                peaks = find(signchanges == -2);
% 
%                if length(peaks) >2
%                    if peaks(1) >7 || peaks(3) < (length(signchanges) - 7)
%                         putativeRippleisARealRipple = 1;
%                    end
%                end
%            end
%         end
         
       if putativeRippleisARealRipple == 1
           realRippleCounter = realRippleCounter+1;
           RealRipples(realRippleCounter) = k;
         
       end
    end
    
    EachRipple = struct();
    
    for z = 1:length(RealRipples)
        Filtered_MedialArrayData_RippleYesNo(i,RippleIndices{RealRipples(z)}) = 1;
        EachRipple(z).Ripple = RippleIndices{RealRipples(z)};
    end    
    MedialArrayRipples(i).EachRipple = EachRipple;
    MedialArrayRipples(i).TimeAxis = TimeAxisM_ds;
    MedialArrayRipples(i).RippleYesNo = Filtered_MedialArrayData_RippleYesNo;
    
    
end

%And lateral array:
Filtered_LateralArrayData_RippleYesNo = zeros(size(LateralArrayData_ds));
LateralArrayRipples = struct;

for i = 1:96

    cutoff = prctile(Filtered_LateralArrayData_rms(i,:),90);     %anywhere the filtered RMS is >90% is putative Ripple
    
    pRipple = find(Filtered_LateralArrayData_rms(i,:)>cutoff);   %These are all the indices of "putative ripple"
    dRipple = diff(pRipple);                                    %this is the number of indices (= milliseconds because 1000 Hz sampling) between Ripple indices
    dRipple = [dRipple 1];
       
    RippleIndices = cell([sum(dRipple>40)+1 1]);
    rippleID = 1;

    
    for j = 1:length(pRipple)
        
        if dRipple(j) > 40                                                      %if the next ripple is more tha 40ms apart, it's a new ripple
        rippleID = rippleID+1;
        end
        RippleIndices{rippleID}  = [RippleIndices{rippleID}, pRipple(j)];       %Add indices onto this Ripple
          
    end
    
    %6) Each putative ripple event was then evaluated based on the number of distinct peaks in the LFP signal (low-passed at 120 Hz) surrounding the event center
    %-- A 40-ms time bin was shifted (5 ms per shift) across ±50 ms, and at least one such time bin must include more than three peaks 
    %-- The first and the last peak cannot be both less than 7 ms away from the edges of this time bin for the putative ripple event to be considered for subsequentanalyses. 
    %-- In addition, the distance between two consecutive ripple centers must exceed 40 ms. 
    %7) To determine the duration of each ripple, eft edge was defined as the first time point before the ripple centerwhere RMS value exceeded the percentile threshold used for initial rip-ple detection, and the right edge was defined as the first time point afterthe ripple center where RMS value fell below the threshold; all othertime points within the duration of this ripple must exceed the RMSthreshold.
    realRippleCounter = 0;
    RealRipples = [];
    
    for k = 1:length(RippleIndices)
        
        putativeRippleindices = RippleIndices{k};
        
        putatitveRippleLP120Hz = LP_Filtered_LateralArrayData_ds(i,putativeRippleindices);
        
        putativeRippleisARealRipple = 0;
        %Need to check for 3 peaks within 40 ms; to have three peaks        
        %a peak will be defined as where the diff of the sign of the diff
        %of the waveform is -2
        %Additionally, if either the first and last peak need to be more than
        %7 ms from the edge of the time bin, then the ripple must be at least 10 ms long 

        %Only consider segments that are at least 41 ms long:
        
        if length(putatitveRippleLP120Hz) > 40        
            for m = 1:length(putatitveRippleLP120Hz)-41
               dPR = diff(putatitveRippleLP120Hz(m:m+40));
               signchanges = diff(sign(dPR));
               peaks = find(signchanges == -2);     %signchanges == -2 are the indices where a peak occurs (where the first differential goes from positive to negative)
               if length(peaks) >2
                   if peaks(1) >7 || peaks(3) < (length(signchanges) - 7)
                        putativeRippleisARealRipple = 1;
                   end
               end
            end
        end    
          
        %%%if I wanted to also consider shorter ripples I could include this
        %code:
%         else
%            if length(putatitveRippleLP120Hz) > 10      
%                dPR = diff(putatitveRippleLP120Hz);            
%                signchanges = diff(sign(dPR));
%                peaks = find(signchanges == -2);
% 
%                if length(peaks) >2
%                    if peaks(1) >7 || peaks(3) < (length(signchanges) - 7)
%                         putativeRippleisARealRipple = 1;
%                    end
%                end
%            end
%         end
         
       if putativeRippleisARealRipple == 1
           realRippleCounter = realRippleCounter+1;
           RealRipples(realRippleCounter) = k;
         
       end
    end
    
    EachRipple = struct();
    
    for z = 1:length(RealRipples)
        Filtered_LateralArrayData_RippleYesNo(i,RippleIndices{RealRipples(z)}) = 1;
        EachRipple(z).Ripple = RippleIndices{RealRipples(z)};
    end    
    LateralArrayRipples(i).EachRipple = EachRipple;
    LateralArrayRipples(i).TimeAxis = TimeAxisL_ds;
    LateralArrayRipples(i).RippleYesNo = Filtered_LateralArrayData_RippleYesNo;
end

fname = ['Ripples/RippleStruct' num2str(exampleScriptnumber) '.mat']
save(fname,'MedialArrayRipples','LateralArrayRipples','-v7.3')


end

end