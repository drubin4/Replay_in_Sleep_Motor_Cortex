%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021

%%This script is going to look overnight at Ripples, and look
%%for a relationship to replay events. 

%This will also incorporate the NS5 data that was separately compiled

%% Open intracranial neural data as well as surface EEG data:

clear all
close all
clc

%% Load the NS5 data

rippledir = 'Ripples';

rippleNS5files = what(rippledir);
rippleNS5files = rippleNS5files.mat;


for i = 1:length(rippleNS5files)
    clear A
    
    putativefilename = rippleNS5files{i};   
    if putativefilename(end-4) ~= 'd'
    tic
    A = load([rippledir '/' putativefilename]);
        
    MedialTime = A.MedialArrayRipples(96).TimeAxis;
    MedialRipple = A.MedialArrayRipples(96).RippleYesNo;
    
    LateralTime = A.LateralArrayRipples(96).TimeAxis;
    LateralRipple = A.LateralArrayRipples(96).RippleYesNo;
    
    putativefilename(end-3:end) = [];
    save([rippledir '/' putativefilename '_simplified.mat'],'MedialTime','MedialRipple','LateralTime','LateralRipple','-v7.3')
    toc
    end
    i
    length(rippleNS5files)
end

