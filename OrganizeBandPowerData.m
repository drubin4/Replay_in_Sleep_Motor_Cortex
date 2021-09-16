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
tic
load('ConcatenatedBandPowerData.mat')
toc
%% Sort the big time axis in order of ascending time

[TimeAxis TimeIndices] = sort(BigTimeAxis,'ascend');

%% Break the data from each band into each array: Array1
tic
Band1_array1 = AllBand1Channels(:,1:96);
Band1_array2 = AllBand1Channels(:,97:192);
Band1_Both_Arrays = AllBand1Channels(:,1:192);

clear AllBand1Channels

Band1_array1 = Band1_array1(TimeIndices,:);
Band1_array2 = Band1_array2(TimeIndices,:);
Band1_Both_Arrays = Band1_Both_Arrays(TimeIndices,:);

save('Band1_Array_1.mat','Band1_array1','TimeAxis','-v7.3')
save('Band1_Array_2.mat','Band1_array2','TimeAxis','-v7.3')
save('Band1_Both_Arrays.mat','Band1_Both_Arrays','TimeAxis','-v7.3')

clear Band1_array1 Band1_array2 Band1_Both_Arrays
toc
%% Array2
tic
Band2_array1 = AllBand2Channels(:,1:96);
Band2_array2 = AllBand2Channels(:,97:192);
Band2_Both_Arrays = AllBand2Channels(:,1:192);

clear AllBand2Channels

Band2_array1 = Band2_array1(TimeIndices,:);
Band2_array2 = Band2_array2(TimeIndices,:);
Band2_Both_Arrays = Band2_Both_Arrays(TimeIndices,:);

save('Band2_Array_1.mat','Band2_array1','TimeAxis','-v7.3')
save('Band2_Array_2.mat','Band2_array2','TimeAxis','-v7.3')
save('Band2_Both_Arrays.mat','Band2_Both_Arrays','TimeAxis','-v7.3')

clear Band2_array1 Band2_array2 Band2_Both_Arrays
toc
%%% Array3
tic
Band3_array1 = AllBand3Channels(:,1:96);
Band3_array2 = AllBand3Channels(:,97:192);
Band3_Both_Arrays = AllBand3Channels(:,1:192);

clear AllBand3Channels

Band3_array1 = Band3_array1(TimeIndices,:);
Band3_array2 = Band3_array2(TimeIndices,:);
Band3_Both_Arrays = Band3_Both_Arrays(TimeIndices,:);

save('Band3_Array_1.mat','Band3_array1','TimeAxis','-v7.3')
save('Band3_Array_2.mat','Band3_array2','TimeAxis','-v7.3')
save('Band3_Both_Arrays.mat','Band3_Both_Arrays','TimeAxis','-v7.3')

clear Band3_array1 Band3_array2 Band3_Both_Arrays
toc
%%% Array4
tic
Band4_array1 = AllBand4Channels(:,1:96);
Band4_array2 = AllBand4Channels(:,97:192);
Band4_Both_Arrays = AllBand4Channels(:,1:192);

clear AllBand4Channels

Band4_array1 = Band4_array1(TimeIndices,:);
Band4_array2 = Band4_array2(TimeIndices,:);
Band4_Both_Arrays = Band4_Both_Arrays(TimeIndices,:);

save('Band4_Array_1.mat','Band4_array1','TimeAxis','-v7.3')
save('Band4_Array_2.mat','Band4_array2','TimeAxis','-v7.3')
save('Band4_Both_Arrays.mat','Band4_Both_Arrays','TimeAxis','-v7.3')

clear Band4_array1 Band4_array2 Band4_Both_Arrays
toc
%%% Array5
tic
Band5_array1 = AllBand5Channels(:,1:96);
Band5_array2 = AllBand5Channels(:,97:192);
Band5_Both_Arrays = AllBand5Channels(:,1:192);

clear AllBand1Channels

Band5_array1 = Band5_array1(TimeIndices,:);
Band5_array2 = Band5_array2(TimeIndices,:);
Band5_Both_Arrays = Band5_Both_Arrays(TimeIndices,:);

save('Band5_Array_1.mat','Band5_array1','TimeAxis','-v7.3')
save('Band5_Array_2.mat','Band5_array2','TimeAxis','-v7.3')
save('Band5_Both_Arrays.mat','Band5_Both_Arrays','TimeAxis','-v7.3')

clear Band5_array1 Band5_array2 Band5_Both_Arrays
toc
