function [cleanECG,kalmanFilterParameters,status,KALMAN_RESET]=kalmanFilter(cleanECG,kalmanFilterParameters,filteredData,peakLocation,instantaneousRR,transformationMatrix,vcg2ecgMatrix, channelsToUse,KALMAN_RESET)


%The Kalman filter filters the filteredData by averaging several
%consecutive ECG complexes. The Kalman gain is updated based on measurement
%and process noise covariances (according to paper by Vullings et al.,
%2011).

%Inputs:
%<cleanECG>: previously filtered ECG complex
%<parameters>: settings, but also variances determined in previous iteration
%<filteredData>: raw data (the name filtered... indicates that it has been bandpass filtered
%<peakLocation>: location of the R-peak in filteredData
%<instantaneousRR>: length between peakLocation and previous one. Used to define length of ECG complex

%Outputs:
%updated parameters and filtered version of current ECG complex

startComplex=peakLocation-round(kalmanFilterParameters.complexStart*kalmanFilterParameters.complexElongation*instantaneousRR)+1;
endComplex=peakLocation-round(kalmanFilterParameters.complexStart*kalmanFilterParameters.complexElongation*instantaneousRR)+round(kalmanFilterParameters.complexElongation*instantaneousRR);
if size(filteredData,2)<endComplex
    %Not enough data available. Go back to main and wait
    status='notFiltered';
        
elseif startComplex<=0
    %do a reset of the filter (and algo)
    KALMAN_RESET=1;
    disp('Kalman bug ==> RESET')
    status='notFiltered';
    
elseif cleanECG.data==zeros(size(cleanECG.data)) & strcmp(cleanECG.firstIteration,'no')==1  
    KALMAN_RESET=1;
    disp('Kalman bug 2 ==> RESET')
    status='notFiltered';
else
    %Define new ECG complex
    newECG=filteredData(:,startComplex:endComplex)'; %transpose as this makes calculations in this function more easy to implement
    if strcmp(cleanECG.firstIteration,'yes')==1        
        filteredECG=zeros(size(newECG));
    else
                    
        firstSampleOfPreviousComplex=cleanECG.lastPeakLocation-round(kalmanFilterParameters.complexStart*kalmanFilterParameters.complexElongation*instantaneousRR)+1;
        lastSampleOfPreviousComplex=cleanECG.lastPeakLocation-round(kalmanFilterParameters.complexStart*kalmanFilterParameters.complexElongation*instantaneousRR)+round(kalmanFilterParameters.complexElongation*instantaneousRR);
        
        %Zero-pad in case cleanECG.data is not long enough (e.g. because of missed peak, current ECG can be twice as long as usual)                        
        filteredECG(1:-firstSampleOfPreviousComplex+1,channelsToUse)=0;
        % AS: 
        %filteredECG(max([-firstSampleOfPreviousComplex+2,1]):max([-firstSampleOfPreviousComplex+2,1])+min([lastSampleOfPreviousComplex,size(cleanECG.data,2)])-max([1,firstSampleOfPreviousComplex]),:)=cleanECG.data(channelsToUse,max([1,firstSampleOfPreviousComplex]):min([lastSampleOfPreviousComplex,size(cleanECG.data,2)]))'; %last ECG complex of the already cleaned ECG. Somewhere (in cleanECG?), the location of the last peak has to be stored
        filteredECG=[filteredECG(:,channelsToUse); cleanECG.data(channelsToUse,max([1,firstSampleOfPreviousComplex]):min([lastSampleOfPreviousComplex,size(cleanECG.data,2)]))'];
        filteredECG(size(filteredECG,1)+1:lastSampleOfPreviousComplex,:)=0; %AS: remove this line?
    end
    newECG(size(filteredECG,1)+1:end,:)=[];
    
    
    kalmanFilterParameters.Sigma=circshift(kalmanFilterParameters.Sigma,[0 -1]);
    
    %Estimate measurement noise (if not specified)
    if strcmp(kalmanFilterParameters.modeMeasurementNoise,'auto')==1
        %Find measurement noise based on non-consistent spatial information
        if strcmp(kalmanFilterParameters.methodMeasurementNoise,'vcg')==1
            %Use VCG to calculate error signal (error signal is part of data
            %that is not mutual between various signal)
            
            errorSignal=newECG-(vcg2ecgMatrix*(transformationMatrix*newECG'))';
            kalmanFilterParameters.Sigma(channelsToUse,kalmanFilterParameters.N)=var(errorSignal,0,1);
            
        elseif strcmp(kalmanFilterParameters.methodMeasurementNoise,'optimal')==1 %optimal in the sense of an optimal estimate of the ECG signal. Not optimal in the sense of the optimal error signal (this approach can underestimate the error signal somewhat when the noise is not white)
            for ch=1:size(newECG,2)
                %Estimate each newECG channel as a linear combination of the
                %other channels: y=x*a y:[Tx1], a:[Nx1], x:[TxN]
                x=newECG(:,[1:ch-1,ch+1:size(newECG,2)]);
                a=(x'*x)\x'*newECG(:,ch);
                estimateECG=x*a;
                errorSignal=newECG(:,ch)-estimateECG;
                kalmanFilterParameters.Sigma(ch,kalmanFilterParameters.N)=var(errorSignal);
            end
        end
    elseif strcmp(kalmanFilterParameters.modeMeasurementNoise,'fixed')==1
        %Fixed measurement noise covariance. This can be based on offline
        %analysis of the data-acquisition system or some other measure.
        
        
        %When set to fixed, the value for <kalmanFilterParameters.Sigma>
        %can have been set by the user in the main.m (or typed into the
        %LabView interface). Anyhow, it is known when this routine is
        %called, so no need to do anything here anymore.
        
    end
    
    
    %Calculate Kalman filter parameters only for actual data, not for
    %zero-padded parts.
    
    if strcmp(cleanECG.firstIteration,'yes')==1
        nonZeros=1:size(filteredECG,1);
    else
        nonZeros=find(sum(filteredECG.^2,2)~=0);
    end
    
    %Calculate rho
    kalmanFilterParameters.rho=circshift(kalmanFilterParameters.rho,[-1 0]);
    kalmanFilterParameters.rho(kalmanFilterParameters.N,channelsToUse)=mean((newECG(nonZeros,:)-filteredECG(nonZeros,:)).^2,1);
    
    %find non-empty heartbeats in rho, meaning that they have been filled in
    %previous iterations
    nonEmptyBeats=find(sum(kalmanFilterParameters.rho.^2,2)~=0);
    lambda=mean(kalmanFilterParameters.rho(nonEmptyBeats,channelsToUse),1)'-mean(kalmanFilterParameters.Sigma(channelsToUse,nonEmptyBeats),2)-kalmanFilterParameters.psi(channelsToUse);
    %lambda should have dimensions [Tx1] or [1xT] with T the number of channels
    %(=8, for ECG filtering, just use all channels)
    lambda(lambda<0)=0; %set to minimally zero
    
    %Calculate Kalman gain
    if strcmp(cleanECG.firstIteration,'yes')==1
        K=ones(size(lambda));
    else
        K=(lambda+kalmanFilterParameters.psi(channelsToUse))./(lambda+kalmanFilterParameters.psi(channelsToUse)+mean(kalmanFilterParameters.Sigma(channelsToUse,nonEmptyBeats),2));       
    end
    
    %Ensure that all channels are filtered equally much. For now, use the
    %mean of the Kalman gain for every channel. Probably a more elegant
    %approach would somehow weigh the channels but this would require a
    %completely new derivation of the Kalman filter equations.
    K=mean(K)*ones(size(K));
    
    
    %For debugging. K seems always near 1, so hardly any filter effect. Has
    %to do with estimation of measurement noise. Might be better if this is
    %fixed (based on location and circumstances of recording)
        %fprintf('K:\n'); disp(K);fprintf('\n\n');
    
    kalmanFilterParameters.psi(channelsToUse)=kalmanFilterParameters.psi(channelsToUse)+lambda-K.*(kalmanFilterParameters.psi(channelsToUse)+lambda);
    
    filteredECG(nonZeros,:)=filteredECG(nonZeros,:)+(newECG(nonZeros,:)-filteredECG(nonZeros,:))*diag(K); %keep zeros from zero-padding. Otherwise for some parts of filtered signal, the newly arriving data is not filtered. Better to show nothing, than suddenly some unfiltered data
    
    status='filtered';
    
    %filteredECG is the filtered version of the incoming ECG complex. It is,
    %however, not a continuous filtered signal. Hence, filteredECG must be
    %appended to the variable cleanECG
    
    %Make sure not to overwrite old data and do not append zeros trailing
    %the filteredECG
    
    
    
    firstSampleFilteredECG=max([cleanECG.lastWrittenSample-startComplex+1,1]);
    lastSampleFilteredECG=find(sum(filteredECG.^2,2)~=0,1,'last');
        
    if isempty(cleanECG.lastWrittenSample)==1
        cleanECG.lastWrittenSample=startComplex; %AS:-1
    end
    %AS: this line creates the shift:
    %cleanECG.data(channelsToUse,cleanECG.lastWrittenSample+1:cleanECG.lastWrittenSample+lastSampleFilteredECG-firstSampleFilteredECG)=filteredECG(firstSampleFilteredECG+1:lastSampleFilteredECG,:)'; %AS: remove the +1
    
    
    %What part of data is overlapping
    newOverlappingSamples=find(sum(filteredECG(1:firstSampleFilteredECG,:).^2,2)~=0); %in filteredECG 
    oldOverlappingSamples=cleanECG.lastWrittenSample-firstSampleFilteredECG+newOverlappingSamples; %in cleanECG %AS: why doing +newOverlappingSamples?
    
    %AS: added by Aline:
    cleanECG.data(channelsToUse,startComplex+numel(newOverlappingSamples):startComplex+numel(newOverlappingSamples)+lastSampleFilteredECG-firstSampleFilteredECG-1)=filteredECG(firstSampleFilteredECG+1:lastSampleFilteredECG,:)';

    
    %Update parameters
    cleanECG.lastPeakLocation=peakLocation;
    %AS:
    %cleanECG.lastWrittenSample=cleanECG.lastWrittenSample+1+lastSampleFilteredECG-firstSampleFilteredECG;
    cleanECG.lastWrittenSample=startComplex+numel(newOverlappingSamples)+lastSampleFilteredECG-firstSampleFilteredECG-1;
    
    %For overlapping data, fade out data from previous heartbeat, fade in data from current heartbeat
    unfadingGradient=1./(1+exp(-kalmanFilterParameters.fadingGradientAlpha*((1:length(oldOverlappingSamples))-length(oldOverlappingSamples)/2)));
    fadingGradient=ones(size(unfadingGradient))-unfadingGradient;
    overlappingData=cleanECG.data(channelsToUse,oldOverlappingSamples).*repmat(fadingGradient,[size(channelsToUse,2) 1])+filteredECG(newOverlappingSamples,:)'.*repmat(unfadingGradient,[size(channelsToUse,2) 1]); %CHANGED !!! cleanECG.data
    
    cleanECG.data(channelsToUse,oldOverlappingSamples)=overlappingData;
    
    %Lowpass filter for transitions
    %for h=1:size(cleanECG.data,1)
    for h=channelsToUse
        cleanECG.data(h,:)=filtfilt(kalmanFilterParameters.filter(1,:),kalmanFilterParameters.filter(2,:),cleanECG.data(h,:));
    end

    cleanECG.firstIteration='no';
end











