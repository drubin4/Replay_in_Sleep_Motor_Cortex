%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021


%% Open intracranial neural data as well as surface EEG data:

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
    

%% First open an intracranial data set of interest:

%DataSetToOpen = 'ncTX_Array_1.mat';
DataSetToOpen = 'ncTX_Array_2.mat';
%DataSetToOpen = 'SpikePower_Array_1.mat';
%DataSetToOpen = 'SpikePower_Array_2.mat';

load(DataSetToOpen)

%SP = single(ncTX_array1);
%clear ncTX_array1
SP = single(ncTX_array2);
clear ncTX_array2
%SP = SpikePower_array1;
%clear SpikePower_array1
%SP = SpikePower_array2;
%clear SpikePower_array2

%% Next open the EEG EDF File:

EDFFilename = '*EDFFILE*.edf';
filename = EDFFilename;

EpilogLoadWithoutDisplay

%To convert time axes from DateNums to DateTimes use "datetime(datevec())"

RawEEGData_DateTimeAxis = datetime(datevec(RawEEGData_TimeAxis));
EEG_Spectrogram_DateTimeAxis = datetime(datevec(EEG_Spectrogram_TimeAxis));

%To convert a DateTime axis to DateNum use "datenum()":

EEGTimeAxis = datenum(RawEEGData_DateTimeAxis);
EEGSpectrogramTimeAxis = datenum(EEG_Spectrogram_DateTimeAxis);
%%
T_start = datetime(2020,11,23,18,0,0);
T_end = datetime(2020,11,25,12,0,0);

%% Plot The Spectrogram
figure(3)
h = imagesc(EEGSpectrogramTimeAxis, EEG_Spectrogram_Frequency_Axis, EEG_Spectrogram);
xlabel('Time (hr:min)');

c_lim = [-20 10];
caxis(c_lim);
f_lim = [0 20];

t_lim = [datenum(T_start) datenum(T_end)];
%t_lim = [-inf inf]

axis tight;
axis([t_lim f_lim]);
datetick('x','keeplimits');

set(gca, 'YDir', 'normal');

ylabel('Freq (Hz)');
xlabel('Time');
title('EEG Spectrogram')

%% Look at the raw EEG and the intracranial data at the time of the first Sync event

% Plot raw EEG of first sync event:

figure(102)
plot(RawEEGData_DateTimeAxis,RawEEGData)
xlim([datetime(2020,11,23,22,58,0) datetime(2020,11,23,23,0,0)])
xlabel('Time')
ylabel('uV')
title('Sync Event #1')


figure(103)
plot(RawEEGData_DateTimeAxis,RawEEGData)
xlim([datetime(2020,11,23,22,59,35) datetime(2020,11,23,22,59,55)])
xlabel('Time')
ylabel('uV')
title('Sync Event #1')

%% %Open the NS5 file from the first sync event:

NS5dir = '*DataDirectory*\Data\_Lateral\NSP Data';
NS5Filename = [NS5dir '\NSP_LATERAL_2020_*****(2)002.ns5'];

NSData_sync = openNSx(NS5Filename);

%Time Start of NS5 Data file
NS5SyncDateTime = NSData_sync.MetaTags.DateTimeRaw;
NS5SyncDateTime(3) = [];            %discard day of week
NS5SyncDateTime(6) = NS5SyncDateTime(6) + NS5SyncDateTime(7)/1000;
NS5SyncDateTime(7) = [];
NS5SyncDateTime = datetime(NS5SyncDateTime);
NS5SyncDateTime = NS5SyncDateTime - hours(5);

%Time Axis of NS5DataFile:

Duration_of_NS5_file =  length(NSData_sync.Data{1})*(1/30000) + length(NSData_sync.Data{2})*(1/30000);
NS5_time_axis = NS5SyncDateTime:seconds((1/30000)):((NS5SyncDateTime+seconds(Duration_of_NS5_file))-seconds(1/30000));

%NS5 sync event axis:
SyncChannel = [NSData_sync.Data{1} NSData_sync.Data{2}];
SyncChannel = SyncChannel(97,:);

%Get rid of first 8 seconds:
SyncChannel(1:240000) = [];
NS5_time_axis(1:240000) = [];

%%
TimeOffset = -0.95;

figure(103)
clf
figure(103)
set(103,'Position',[50 100 1000 600])

subplot(2,2,1)
hold on
plot(NS5_time_axis,SyncChannel/15,'r')
plot(RawEEGData_DateTimeAxis,RawEEGData,'b')
xlim([datetime(2020,11,23,22,58,30) datetime(2020,11,23,22,59,0)])
ylim([-300 300])
xlabel('Time')
ylabel('uV')
title('Before Adjusting')


subplot(2,2,2)
hold on
plot(NS5_time_axis,SyncChannel/20,'r')
plot(RawEEGData_DateTimeAxis,RawEEGData,'b')
xlim([datetime(2020,11,23,22,58,35) datetime(2020,11,23,22,58,40)])
ylim([-300 300])
xlabel('Time')
ylabel('uV')
title('Before Adjusting')

subplot(2,2,3)
hold on
plot(NS5_time_axis,SyncChannel/20,'r')
plot(RawEEGData_DateTimeAxis-seconds(TimeOffset),RawEEGData,'b')
xlim([datetime(2020,11,23,22,58,30) datetime(2020,11,23,22,59,0)])
ylim([-300 300])
xlabel('Time')
ylabel('uV')
title('After Adjusting')

subplot(2,2,4)
hold on
plot(NS5_time_axis,SyncChannel/15,'r')
plot(RawEEGData_DateTimeAxis-seconds(TimeOffset),RawEEGData,'b')
xlim([datetime(2020,11,23,22,58,35) datetime(2020,11,23,22,58,40)])
ylim([-300 300])
xlabel('Time')
ylabel('uV')
title('After Adjusting')

Title_Details = sprintf(', Offset Time: %4.0f ms',TimeOffset*1000);
suptitle(['Sync Event #1' Title_Details])


fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 18 11];
fig.PaperSize = [18 11];
%print('Syncing_EEG_and_Intracranial_Data.pdf','-bestfit','-dpdf')



%% Plot the second sync event
%% Look at the raw EEG and the intracranial data at the time of the second Sync event

%Second sync occurred at: block 66) 09:44:03 

% Plot raw EEG of first sync event:

figure(102)
plot(RawEEGData_DateTimeAxis,RawEEGData)
xlim([datetime(2020,11,24,9,43,45) datetime(2020,11,24,9,44,30)])
xlabel('Time')
ylabel('uV')
title('Sync Event #2')


figure(103)
plot(RawEEGData_DateTimeAxis,RawEEGData)
xlim([datetime(2020,11,24,9,44,0) datetime(2020,11,24,9,44,5)])
xlabel('Time')
ylabel('uV')
title('Sync Event #2')

%% Open the NS5 file from the second sync event:

NS5Filename = [NS5dir '\NSP_LATERAL_2020_1124_******3(66)129.ns5'];

NSData_sync = openNSx(NS5Filename);

%Time Start of NS5 Data file
NS5SyncDateTime = NSData_sync.MetaTags.DateTimeRaw;
NS5SyncDateTime(3) = [];            %discard day of week
NS5SyncDateTime(6) = NS5SyncDateTime(6) + NS5SyncDateTime(7)/1000;
NS5SyncDateTime(7) = [];
NS5SyncDateTime = datetime(NS5SyncDateTime);
NS5SyncDateTime = NS5SyncDateTime - hours(5);

%Time Axis of NS5DataFile:

Duration_of_NS5_file =  length(NSData_sync.Data{1})*(1/30000) + length(NSData_sync.Data{2})*(1/30000);
NS5_time_axis = NS5SyncDateTime:seconds((1/30000)):((NS5SyncDateTime+seconds(Duration_of_NS5_file))-seconds(1/30000));

%NS5 sync event axis:
SyncChannel = [NSData_sync.Data{1} NSData_sync.Data{2}];
SyncChannel = SyncChannel(97,:);

%Get rid of first 8 seconds which has other mark in it:
SyncChannel(1:240000) = [];
NS5_time_axis(1:240000) = [];

%%

TimeOffset = -2.5;

figure(1003)
clf
figure(1003)
set(1003,'Position',[50 100 1000 600])

subplot(2,2,1)
hold on
plot(NS5_time_axis,SyncChannel/15,'r')
plot(RawEEGData_DateTimeAxis,RawEEGData,'b')
xlim([datetime(2020,11,24,9,44,30) datetime(2020,11,24,9,45,0)])
ylim([-1000 1000])
xlabel('Time')
ylabel('uV')
title('Before Adjusting')


subplot(2,2,2)
hold on
plot(NS5_time_axis,SyncChannel/15,'r')
plot(RawEEGData_DateTimeAxis,RawEEGData,'b')
xlim([datetime(2020,11,24,9,44,36) datetime(2020,11,24,9,44,44)])
ylim([-1000 1000])
xlabel('Time')
ylabel('uV')
title('Before Adjusting')

subplot(2,2,3)
hold on
plot(NS5_time_axis,SyncChannel/15,'r')
plot(RawEEGData_DateTimeAxis-seconds(TimeOffset),RawEEGData,'b')
xlim([datetime(2020,11,24,9,44,30) datetime(2020,11,24,9,45,0)])
ylim([-1000 1000])
xlabel('Time')
ylabel('uV')
title('After Adjusting')

subplot(2,2,4)
hold on
plot(NS5_time_axis,SyncChannel/15,'r')
plot(RawEEGData_DateTimeAxis-seconds(TimeOffset),RawEEGData,'b')
xlim([datetime(2020,11,24,9,44,36) datetime(2020,11,24,9,44,44)])
ylim([-1000 1000])
xlabel('Time')
ylabel('uV')
title('After Adjusting')


Title_Details = sprintf(', Offset Time: %4.0f ms',TimeOffset*1000);
suptitle(['Sync Event #2' Title_Details])


fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 18 11];
fig.PaperSize = [18 11];
print('Syncing_EEG_and_Intracranial_Data_2.pdf','-bestfit','-dpdf')


%% Sync time axes and then save EEG data with adjusted Time Axes:

RawEEGData_DateTimeAxis = RawEEGData_DateTimeAxis - seconds(TimeOffset);
EEG_Spectrogram_DateTimeAxis = EEG_Spectrogram_DateTimeAxis - seconds(TimeOffset);

save('Synced_EEG_Data.mat','RawEEGData','EEG_Spectrogram','EEG_Spectrogram_Frequency_Axis','RawEEGData_DateTimeAxis','EEG_Spectrogram_DateTimeAxis')


