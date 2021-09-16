%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021


%% Concatenate all SpikePower and ncTX spike trains from separate SLC files and read NS5 files to get absolute time for comparison with EEG data

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

SLCdir = '*DefineDataDirectoryHere*\Data\SLCData';
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
        
    %Record the neural data from the SLC file    
    
    if i == 1
        AllSpikePowerChannels = SLCFile.spikePower.values;
        AllncTXChannels = SLCFile.ncTX.values;
        
    end
    if i > 1           
        AllSpikePowerChannels = [AllSpikePowerChannels; SLCFile.spikePower.values];
        AllncTXChannels = [AllncTXChannels; SLCFile.ncTX.values];
    end
    
    
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
       
           
               if strcmp(tempFname(13:end-7),unique_portion_of_ID) && strcmp(tempFname(end-2:end),'ns3')

                   DateTimeMatch = 'MATCH'
                   SLC_filenames{i}
                   NS5_filename = tempFname
               end
           
           
           %Need to create a special rule for the overnight files that have a
           %different naming convention

%            if strcmp(unique_portion_of_ID(1:19),'2020_0225_000955(5)')
% 
%                Overnight_file_code = SLC_filenames{i};
% 
%                %pull the 3 digit number from the end of the filename
%                Overnight_file_code = Overnight_file_code(end-7:end-5);
%                %convert "f" to "0" to match NS3/NS5 files
%                if strcmp(Overnight_file_code(1),'f')
%                    Overnight_file_code(1) = '0';
%                end
% 
%                if strcmp(tempFname(1:15),'20200225-000730') && strcmp(tempFname(end-6:end-4),Overnight_file_code)
% 
%                    DateTimeMatch = 'OVERNIGHT MATCH'
%                    SLC_filenames{i}
%                    NS5_filename = tempFname
% 
%                end
% 
%            end

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
        
        %timing is stored differently for the overnight sessions:
%         if strcmp(DateTimeMatch,'OVERNIGHT MATCH')
%             %Use the nsp2Clock from the SLC file as the Time axis:
%                 TempTimeAxisSeconds = SLCFile.clocks.nsp2Clock;
% 
%                     %Look for discontinuity in time axis from NSP sync process:
%                     diff_TTA = diff(TempTimeAxisSeconds);
%                     resettime = find(diff_TTA<0);
%                     if ~isempty(resettime)
%                         for jj = 1:length(resettime)
%                             TempTimeAxisSeconds((resettime(jj)+1):end) = TempTimeAxisSeconds((resettime(jj)+1):end) + TempTimeAxisSeconds(resettime(jj));               
%                         end
%                     end
%                     TempTimeAxis = StartDateTime + seconds(TempTimeAxisSeconds);
%                     TempTimeAxis = TempTimeAxis';
%             
%             
%         end
            
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
save('ConcatenatedData.mat','AllSpikePowerChannels','AllncTXChannels','BigTimeAxis','-v7.3')
    