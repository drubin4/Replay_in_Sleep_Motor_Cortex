%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021


function [h] = plotEEGSpectrogram(t_start,t_end,figurenumber,EEGSpectrogramTimeAxis,EEG_Spectrogram_Frequency_Axis,EEG_Spectrogram)
%plotEEGSpectrogram: Plot the EEG spectrogram into defined time paramters
%and figure number
EEG_Spectrogram_DateTimeAxis = datetime(datevec(EEGSpectrogramTimeAxis));

StartTimeIndex_EEG = find(EEG_Spectrogram_DateTimeAxis > t_start);
StartTimeIndex_EEG = StartTimeIndex_EEG(1);

EndTimeIndex_EEG = find(EEG_Spectrogram_DateTimeAxis > t_end);
EndTimeIndex_EEG = EndTimeIndex_EEG(1);

figure(figurenumber)
h = imagesc(EEGSpectrogramTimeAxis(StartTimeIndex_EEG:EndTimeIndex_EEG), EEG_Spectrogram_Frequency_Axis, EEG_Spectrogram(:,(StartTimeIndex_EEG:EndTimeIndex_EEG)));
%xlabel('Time (hr:min)');

c_lim = [-20 10];
caxis(c_lim);
f_lim = [0 20];

t_lim = [datenum(t_start) datenum(t_end)];
%t_lim = [-inf inf]

axis tight;
axis([t_lim f_lim]);
datetick('x','keeplimits');
set(gca, 'YDir', 'normal');

ylabel('Freq (Hz)');
%xlabel('Time');
%title('EEG Spectrogram')


end

