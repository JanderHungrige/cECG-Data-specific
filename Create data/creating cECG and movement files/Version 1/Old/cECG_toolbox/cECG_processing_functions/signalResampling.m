function [data,extraSamples,numberNewSamples]=signalResampling(rawData,data,filterCoefficients,numberNewSamples,extraSamples,FS_RATIO)

%Note: data WILL BE DELAYED BY 3/250=12ms

overlap=filterCoefficients.resamplingDelay2; %3 samples overlap at 250Hz
extraSamplesOld=extraSamples; %is counted as a new sample
extraSamples=mod(numberNewSamples+extraSamplesOld,FS_RATIO);
resampledData=resampleData(rawData(:,end-numberNewSamples-extraSamplesOld-filterCoefficients.resamplingDelay2-overlap+1:end-extraSamples),filterCoefficients,FS_RATIO); %take also the delay of previous calls!!
numberNewSamples=(numberNewSamples+extraSamplesOld-extraSamples)/FS_RATIO;
%data=[data(:,numberNewSamples+1:end),resampledData(:,filterCoefficients.resamplingDelay+1:end)]; %compensate for overlap
data=[data(:,numberNewSamples+1:end),zeros(size(data,1),numberNewSamples)]; 
data(:,end-numberNewSamples+1:end)=resampledData(:,overlap/FS_RATIO+1:end); 

end


function resampledData=resampleData(data,filterCoefficients,FS_RATIO)

% data must be a multiple of 32 for resampling!!!! => delay as well
% Resample data from 8kHz to 250Hz

%Low-pass data
filteredData=filter(filterCoefficients.resampling,1,data,[],2);

%Omit samples
resampledData=filteredData(:,1:FS_RATIO:end);

%Delay of half the filter length
resampledData(:,1:filterCoefficients.resamplingDelay)=[]; 

%here resample data is filtered, downsampled and synchronized with the
%input data BUT the 3 first samples are corrupted --> overlap needed

end
