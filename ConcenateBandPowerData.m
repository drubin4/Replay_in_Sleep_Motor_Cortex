%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021


%% Concatenate all LFP data (in each power band) from separate SLC files and read NS5 files to get absolute time for comparison with EEG data

clear all
close all
clc

%% Add all directories to path to allow for ease of access to SLC and NS5 files
% Make "SecondSession" the default directory

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
    
    
%% Generate a list of ALL the SLC files and store it in cell array "SLC_filenames":
clear SLC_filenames

SLCdir = '*DefineDataDirectoryHere*\SLCData';
S = what(SLCdir);
S = S.mat;

j = 1;
for i = 1:length(S)
    putativefilename = S{i};
    
    if strcmp(putativefilename(1:3),'SLC');% || strcmp(putativefilename(1:3),'Int')
        SLC_filenames{j} = putativefilename;
        j = j+1;
    end
end

%% Open each SLC file: 

%Store the time axis here:
BigTimeAxis = [];

%Store the spike power and ncTX data
AllSpikePowerChannels = [];
AllncTXChannels = [];

%In the future consider going back and get the bandpower

%Keep track of any missing data file times
missingTimes = 0;
specificMissingTimes = [];

% Also open each NS5 file header to get the time of the SLC file start
NS5dir = '*DefineDataDirectoryHere*\Data\_Lateral\NSP Data';
NSFiles = dir(NS5dir);


%Open each SLC File:

for i = 1:length(SLC_filenames)
    
% Load the SLC File:
    fname = [SLCdir '\' SLC_filenames{i}];
    tic
    SLCFile = load(fname);
        
    %Record the bandpower data from the SLC file    
%     
    if i == 1
        AllBand1Channels = SLCFile.STFT.values(:,1:192);
        AllBand2Channels = SLCFile.STFT.values(:,193:384);
        AllBand3Channels = SLCFile.STFT.values(:,385:576);
        AllBand4Channels = SLCFile.STFT.values(:,577:768);
        AllBand5Channels = SLCFile.STFT.values(:,769:960);        
    end
    if i > 1           
        AllBand1Channels = [AllBand1Channels; SLCFile.STFT.values(:,1:192)];
        AllBand2Channels = [AllBand2Channels; SLCFile.STFT.values(:,193:384)];
        AllBand3Channels = [AllBand3Channels; SLCFile.STFT.values(:,385:576)];
        AllBand4Channels = [AllBand4Channels; SLCFile.STFT.values(:,577:768)];
        AllBand5Channels = [AllBand5Channels; SLCFile.STFT.values(:,769:960)]; 
    end
    
    FreqBands = GetSTFTBands(SLCFile.sSLC);
    %partition off the unique portion of the SLC file name to find the
    %associated NS5 file based on the naming convention
    
    %Also need to distinguish between the blocks (which start as "SLC" and the
    %"Interblocks" which start as "Interblock"
    
    unique_portion_of_ID = SLC_filenames{i}; 
    
    if strcmp(SLC_filenames{i}(1:3),'Int')
    unique_portion_of_ID = unique_portion_of_ID(20:end-4);
    blocktype = 'Interblock';
    end
    if strcmp(SLC_filenames{i}(1:3),'SLC')
    unique_portion_of_ID = unique_portion_of_ID(9:end-4);    
    blocktype = 'Block';
    end
    
    
    NS5_filename = [];
    %find the associated NS5 file: 
    if strcmp('Block',blocktype)
     1
     for j = 1:length(NSFiles)
       tempFname = NSFiles(j).name;  
       
       if length(tempFname) > 2
       
           
               if strcmp(tempFname(13:end-7),unique_portion_of_ID) && strcmp(tempFname(end-2:end),'ns5')

                   DateTimeMatch = 'MATCH'
                   SLC_filenames{i}
                   NS5_filename = tempFname
               end
           
        end
       
    end
    
    %open the NS5 file and find the start time
    if ~isempty(NS5_filename)
        FullNS5Name = [NS5dir '\' NS5_filename];
        NS5Data = openNSx(FullNS5Name);

        %Record Start Time
        StartTime = NS5Data.MetaTags.DateTimeRaw;

        %Adjust Time Zone from UTC to EST:
       
        StartDateTime = datetime(StartTime(1),StartTime(2),StartTime(4),StartTime(5),StartTime(6),StartTime(7),StartTime(8));
        StartDateTime = StartDateTime - hours(5);
        
    else
        
        StartDateTime = BigTimeAxis(end) + milliseconds(20);
        missingTimes = missingTimes + 1;
        specificMissingTimes(missingTimes) = i;
    end
    
    
        %Create a corresponding time axis using the start time from the NS5
        %file and the Clocks vector from the SLC file:
         
        %for daytime sessions do it one way:
        
        if strcmp(DateTimeMatch,'MATCH')
                %Use the nsp2Clock from the SLC file as the Time axis:
                TempTimeAxisSeconds = SLCFile.clocks.nsp2Clock;

                    %Look for discontinuity in time axis from NSP sync process:
                    diff_TTA = diff(TempTimeAxisSeconds);
                    resettime = find(diff_TTA<0);
                    if ~isempty(resettime)
                        TTAS = TempTimeAxisSeconds;
                        for k = 1:length(resettime)                    
                            TTAS((resettime(k)+1):end) = TTAS((resettime(k)+1):end) + TempTimeAxisSeconds(resettime(k));                
                        end
                        TempTimeAxisSeconds = TTAS;
                    end
                    TempTimeAxis = StartDateTime + seconds(TempTimeAxisSeconds);
                    TempTimeAxis = TempTimeAxis';
        end
        

            
        if i == 1
            BigTimeAxis = TempTimeAxis;
        end
        if i > 1   
            BigTimeAxis = [BigTimeAxis TempTimeAxis];    
        end
        
    
    
    
    
    toc
    i
    end
    
end
save('ConcatenatedBandPowerData.mat','AllBand1Channels','AllBand2Channels','AllBand3Channels','AllBand4Channels','AllBand5Channels','FreqBands','BigTimeAxis','-v7.3')
    