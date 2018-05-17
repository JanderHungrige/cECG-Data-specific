function filteredData = bandpassFiltering(filteredData, data, filterCoefficients, numberNewSamples, normal)

if normal==1
    % NB: this filter don't seem to work on the dataset '20120312_170044 adult on new neonat mat 0dB' when NOT taking differential leads

    overlap=filterCoefficients.bandpassDelay;
%     delay=filterCoefficients.highpassDelay+filterCoefficients.lowpassDelay;
%     newFiltData=filter(filterCoefficients.lowpass,1,data(:,end-numberNewSamples-delay-overlap+1:end),[],2);
%     newFiltData=filter(filterCoefficients.highpass,1,newFiltData,[],2);
    delay=filterCoefficients.bandpassDelay;
    newFiltData=filter(filterCoefficients.bandpass,1,data(:,end-numberNewSamples-delay-overlap+1:end),[],2);
else %sharp bandpass (0) for infant data
    overlap=filterCoefficients.sharpBandpassDelay;
    delay=filterCoefficients.sharpBandpassDelay;
    newFiltData=filter(filterCoefficients.sharpBandpass,1,data(:,end-numberNewSamples-delay-overlap+1:end),[],2);
end

filteredData=[filteredData(:,numberNewSamples+1:end),newFiltData(:,overlap+delay+1:end)];
%filteredData=[filteredData(:,numberNewSamples+1:end),zeros(size(filteredData,1),numberNewSamples)];
%filteredData(:,end-numberNewSamples+1:end)=newFiltData(:,overlap+delay+1:end);