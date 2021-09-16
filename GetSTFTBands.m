%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021


function freqMinMax = GetSTFTBands(sSLC)

if isfield(sSLC, 'sSLC')
    sSLC = sSLC.sSLC;
end


for ii = 1:sSLC.features.STFT.bands
    bandName = sprintf('FreqBins%d', ii);
    band = sSLC.features.STFT.(bandName);
    if all(diff(band)~=1)
        warning('STFT bands %d does not contain contiguous frequencies, %s', num2str( band ) );
    end
    freqMinMax(ii,:) = sSLC.features.STFT.frequencies( [min(band) max(band)] );
end

end
