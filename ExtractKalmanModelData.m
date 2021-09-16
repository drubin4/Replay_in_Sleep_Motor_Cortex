%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021


%% Open an SLC file from during the game to extract the parameters used by the model/Kalman filter

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
    
%% Generate a list of ALL the SLC files and store it in cell array "SLC_filenames":
clear SLC_filenames

SLCdir = '*DataDirectory*\Data\SLCData';
S = what(SLCdir);
S = S.mat;

j = 1;
for i = 1:length(S)
    putativefilename = S{i};
    
    if strcmp(putativefilename(1:3),'SLC')  % || strcmp(putativefilename(1:3),'Int')
        SLC_filenames{j} = putativefilename;
        j = j+1;
    end
end


%% Open an SLC file from during the game

%File #80 works well enough

 fname = [SLCdir '\' SLC_filenames{80}];
    tic
    SLCFile = load(fname);

    KalmanModel = SLCFile.sSLC.decoders.kalman;
    
    A = KalmanModel.A;
    H = KalmanModel.H;
    K = KalmanModel.K;
    SpikePowerMask = KalmanModel.spikePowerMask;
    ncTXMax = KalmanModel.ncTXMask;
    
 %%
    
    save('KalmanModelData.mat','KalmanModel','-v7.3')
    