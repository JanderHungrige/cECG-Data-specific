function filteredData = notchFiltering(filteredData, data, filterCoefficients, numberNewSamples)

overlap=filterCoefficients.notchDelay;
delay=filterCoefficients.notchDelay;
newFiltData=filter(filterCoefficients.notch,1,data(:,end-numberNewSamples-delay-overlap+1:end),[],2);

%filteredData=[filteredData(:,numberNewSamples+1:end),newFiltData(:,overlap+delay+1:end)];
filteredData=[filteredData(:,numberNewSamples+1:end),zeros(size(filteredData,1),numberNewSamples)];
filteredData(:,end-numberNewSamples+1:end)=newFiltData(:,overlap+delay+1:end);


