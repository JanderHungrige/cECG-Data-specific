function transformedData=signalTransformation2(projectedData,transformedData,filterCoefficients,numberNewSamples)

% The size of the vectors <projectedData> must be higher or equal to <filterCoefficients.overlapForSigTransform>
% size(projectedData,2) must be strictly higher than <numberNewSamples>
% otherwise all new samples cannot be buffered without lost

%Transform all new samples
waveletOutput1=waveletAnalysis(projectedData(1,end-numberNewSamples-filterCoefficients.lengthWavelet+1:end),filterCoefficients.wavelet); %take all new samples + waveletLength samples in the past to compensate for 1/2 waveletLength inaccurracy of last call and 1/2 waveletLength inaccuracy to come

%Shift <transformedData>
transformedData=[transformedData(numberNewSamples+1:end),zeros(1,numberNewSamples)];
%Updates the last samples of <transformedData>
transformedData(end-numberNewSamples-ceil(filterCoefficients.lengthWavelet/2)+1:end-ceil(filterCoefficients.lengthWavelet/2)+1)=abs(waveletOutput1); %overwrite on waveletLength/2 past samples (to compensate for inaccuracy of previous call of this function)
% NB: there is a delay introduced due to the concolution with the wavelet: the last
% ceil(filterCoefficients.overlapForSigTransform/2) samples have zero
% values. (=inaccuracy of this call of the function, to be compensate
% during the next call)

end


%Perform wavelet analysis
function WaveletOutput=waveletAnalysis(Signal,wavelet)

TempWaveletOutput=conv(Signal(1,:),wavelet);
%WaveletOutput(1,1+floor(size(wavelet,1)/2):length(Signal)-ceil(size(wavelet,1)/2)+1)=TempWaveletOutput(size(wavelet,1):end-size(wavelet,1)+1);
WaveletOutput=TempWaveletOutput(size(wavelet,1):end-size(wavelet,1)+1);
end