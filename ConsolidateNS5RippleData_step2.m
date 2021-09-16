%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021


%%This script will open the "simplified" NS5 data and will attempt to
%%concatenate all the data into a single LONG pair of matrices

%% Open intracranial neural data as well as surface EEG data:

clear all
close all
clc

%% Load the NS5 data

rippledir = 'Ripples';

rippleNS5files = what(rippledir);
rippleNS5files = rippleNS5files.mat;

MedialArrayRippleYesNo = []; 
LateralArrayRippleYesNo = [];
MedialTimeAxis = [];
LateralTimeAxis = [];    
    
for i = 1:length(rippleNS5files)
    
    putativefilename = rippleNS5files{i};   
    if putativefilename(end-4) == 'd'           
    tic
    A = load([rippledir '/' putativefilename]);
        
    MedialTime_i = A.MedialTime;
    MedialRipple_i = A.MedialRipple;
        
    LateralTime_i = A.LateralTime;
    LateralRipple_i = A.LateralRipple;
    
    
    
    MedialArrayRippleYesNo = [MedialArrayRippleYesNo MedialRipple_i];
    LateralArrayRippleYesNo = [LateralArrayRippleYesNo LateralRipple_i];
    MedialTimeAxis = [MedialTimeAxis MedialTime_i];
    LateralTimeAxis = [LateralTimeAxis LateralTime_i];
        
    i   
    toc
    
    
    end
end

%%
save('Ripples/MedialArrayAllRipples.mat','MedialArrayRippleYesNo','MedialTimeAxis','-v7.3')
save('Ripples/LateralArrayAllRipples.mat','LateralArrayRippleYesNo','LateralTimeAxis','-v7.3')

MedialArrayRippleYesNo = uint8(MedialArrayRippleYesNo);
LateralArrayRippleYesNo = uint8(LateralArrayRippleYesNo);

save('Ripples/MedialArrayAllRipples_lite.mat','MedialArrayRippleYesNo','MedialTimeAxis','-v7.3');
save('Ripples/LateralArrayAllRipples_lite.mat','LateralArrayRippleYesNo','LateralTimeAxis','-v7.3');


