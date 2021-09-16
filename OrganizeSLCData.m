%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021


%% Open the concatenated SLC data, and re-order everything in order of time, then save each array data separately

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
    
    
%% Open the concatenated SLC Data file

load('ConcatenatedData.mat')

%% Sort the big time axis in order of ascending time

[TimeAxis TimeIndices] = sort(BigTimeAxis,'ascend');

%% Break the spike power data and ncTX into separate matrices:

SpikePower_array1 = AllSpikePowerChannels(:,1:96);
SpikePower_array2 = AllSpikePowerChannels(:,97:192);

SpikePower_Both_Arrays = AllSpikePowerChannels(:,1:192);

clear AllSpikePowerChannels

ncTX_array1 = AllncTXChannels(:,1:96);
ncTX_array2 = AllncTXChannels(:,97:192);

ncTX_Both_Arrays = AllncTXChannels(:,1:192);

clear AllncTXChannels

%% Re-order all of the neural data in ascending time

ncTX_array1 = ncTX_array1(TimeIndices,:);
ncTX_array2 = ncTX_array2(TimeIndices,:);
ncTX_Both_Arrays = ncTX_Both_Arrays(TimeIndices,:);

SpikePower_array1 = SpikePower_array1(TimeIndices,:);
SpikePower_array2 = SpikePower_array2(TimeIndices,:);
SpikePower_Both_Arrays = SpikePower_Both_Arrays(TimeIndices,:);

%% Save the data:

save('ncTX_Array_1.mat','ncTX_array1','TimeAxis','-v7.3')
save('ncTX_Array_2.mat','ncTX_array2','TimeAxis','-v7.3')
save('SpikePower_Array_1.mat','SpikePower_array1','TimeAxis','-v7.3')
save('SpikePower_Array_2.mat','SpikePower_array2','TimeAxis','-v7.3')

save('ncTX_Both_Arrays.mat','ncTX_Both_Arrays','TimeAxis','-v7.3')
save('SpikePower_Both_Arrays.mat','SpikePower_Both_Arrays','TimeAxis','-v7.3')

