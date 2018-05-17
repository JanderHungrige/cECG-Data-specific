function filterCoefficients=determineFilters(Fs, adult)
filterCoefficients=struct('lowlowpass',[],'lowlowpassDelay',[],'lowpass',[],'lowpassDelay',[], 'bandpass',[],'bandpassDelay',[],'sharpBandpass',[],'sharpBandpassDelay',[],'highpassDelay',[],'highpass',[],'resampling',[],'resamplingDelay',[],'resamplingDelay2',[],'resamplingLength',[],'sensorTransferFunction',[],'wavelet',[],'lengthWavelet',[]);

%Frequency selective filters to suppress baseline-drift, powerline
%interference, high-frequency components (+ powerline harmonics)



%Temporal filter that improves match between capacitive sensors and
%reference (=galvanic) ECG. This filter was designed to be applied on a
%single signal obtained from spatially combining the various recorded cECG
%signals. Due to linearity of both the spatial and temporal filtering, this
%step can already be applied on the incoming data as to select proper
%signals based on data that already resembles the final goal more.
%load('cECG2gECGfilter.mat');
%filterCoefficients.sensorTransferFunction=cECG2gECGfilter; %for now, filter has not been determined yet

% load('resamplingFilter')
% filterCoefficients.resampling=resamplingFilter.coefficients;
% filterCoefficients.resamplingLength=length(resamplingFilter.coefficients);
% filterCoefficients.resamplingDelay=11;
% filterCoefficients.resamplingDelay2=11*32;

load('resamplingFilter_V2') %Low pass equiripple, 192 coef at 8000Hz, 6 at 250Hz, 45, 205, 0.1, 97
filterCoefficients.resampling=coefficients;
filterCoefficients.resamplingLength=length(filterCoefficients.resampling);
filterCoefficients.resamplingDelay=filterCoefficients.resamplingLength/(8000/Fs)/2;  %3
filterCoefficients.resamplingDelay2=filterCoefficients.resamplingLength/2;  %96 (MULTIPLE OF 32!!)

% if ~adult
%     load('lowlowpass3_Fs250') %106 4 10 1 80 %232 4 7 0.5 80
%     filterCoefficients.lowlowpass=coefficients;
%     filterCoefficients.lowlowpassDelay=ceil(length(filterCoefficients.lowlowpass)); % do overlap instead!
%     
%     load('bandpass10183542_Fs250') %131
%     filterCoefficients.sharpBandpass=Numerator;
%     filterCoefficients.sharpBandpassDelay=ceil(length(filterCoefficients.sharpBandpass)/2);
% end

% load('lowpass35_Fs250')
% filterCoefficients.lowpass=lowpass35_Fs250.Numerator;
% filterCoefficients.lowpassDelay=ceil(length(lowpass35_Fs250.Numerator)/2);
% 
% load('highpass3_Fs250.mat')
% filterCoefficients.highpass=highpass3_Fs250_2.Numerator;
% filterCoefficients.highpassDelay=ceil(length(highpass3_Fs250_2.Numerator)/2);

load('bandpass_Fs250') %Band passequiripple, 292 coef at 250Hz, 1.3-3, 35-38, 60,1,80
filterCoefficients.bandpass=coefficients;
filterCoefficients.bandpassDelay=length(filterCoefficients.bandpass)/2;

%load('notch21242528_Fs250_194') 
load('notch23242526_Fs250_580') %86%
%load('notch232424525_Fs250_1064')
%load('notch2323524525_Fs250_1320') %86%
filterCoefficients.notch=Num;
filterCoefficients.notchDelay=ceil(length(filterCoefficients.notch)/2);

%Filters for signal enhancement in the peak detection. This filter will
%completely distort the ECG, but remains the periodicity. The output of
%this filter will be used for detecting the heart rate, but not use for
%further ECG analysis.
waveletAccuracy=4;
Fc=[23.2963, 30];%17.4722];
waveletBasewidths=Fs*0.2516./Fc;
for i=1:length(waveletBasewidths)
    x=-waveletAccuracy:1/waveletBasewidths(i):waveletAccuracy;
    wavelet=((2/sqrt(3)*pi^(-1/4))*(1-x.^2).*exp(-(x.^2)/2))';
    eval(['wavelet',num2str(i),'=wavelet./ones(length(wavelet),1)*sum(abs(wavelet));']);
end
if adult
    filterCoefficients.wavelet=wavelet1; %wavelet2 for infant!
else
    filterCoefficients.wavelet=wavelet2; 
end
filterCoefficients.lengthWavelet=length(filterCoefficients.wavelet);
%filterCoefficients.overlapForSigTransform=max([length(filterCoefficients.peakDetection1),length(filterCoefficients.peakDetection2)]);
