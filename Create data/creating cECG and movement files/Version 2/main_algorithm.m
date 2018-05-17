% MAIN_ALGORITHM.M
%
% This file:
% - processes the data based on a sliding of size WINDOW_SIZE;
% - creates a video of the signals (if VIDEO=1);
% - creates dynamic graphs of the data processing (if GRAPH=1);
% - saves variables and signals for further use (if SAVE=1);
% - extracts the respiration signal from the 1 kHz injected signal (if RESP=1)
% - computes the amount of motion based on the camera recordings (if CAMERA=1)
% - computes the amount of motion based on the 1 kHz signal (AS: To Be Done!)
%
% This file:
% - is called by launch_processing.m for offline data analysis;
% - is called by LabView for the online demo;
% - calls determineFilters.m, findTransformationMatrix.m, initializeFilter.m, template_dynamic_graph.m,
% template_video_resp.m, template_video.m, hilbertdemod.m, channelrangecorr.m, smooth.m, iirFilter.m,
% signalResampling.m, signalTransformation2.m, peakDetect2.m, kalmanFilter.m, MotionEstimator.m (AS: clean them all!)
%
% The processing of the data consists in different steps:
% - channel selection
% - downsampling
% - temporal filter to compensate for the channel imbalance (AS: To Be Implemented!!)
% - bandpass filtering
% - VCG computation
% - defining the body orientation from the VCG (AS: To Be Implemented!!)
% - Einthoven leads projection
% - peak detection (AS: To Be Improved)
% - Kalman filtering based on the detected peaks for ECG signal enhancement (AS: To Be optimized)
%
% aline.serteyn@gmail.com, r.vullings@tue.nl
% June 2012



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         INITIALIZATION                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if KALMAN_RESET==1
    if DO_INIT %Never reset the following variables/parameters
        
        % Parameters
        % ----------
        if ONLINE_DEMO
            THRESHOLD=0.07;
            MIN_THRESHOLD=0; % change it to remove broken sensors
        else
            THRESHOLD=150000; %300000 originally (video)
            MIN_THRESHOLD=0; % change it to remove broken sensors
        end
        NUMBER_BEFORE_RESET=2; %number of consecutive unreliable windows before a reset is generated due to bad coverage
        FS=500; %was 250 originally
        if ~ADULT
            HR=[55,250];
            meanRR2=100;
            NUM_RR_AVERAGED=20; %this is in number of RR intervals
        else
            HR=[32,200];
            meanRR2=200;
            NUM_RR_AVERAGED=10; %this is in number of RR intervals
        end
        
        if ~FIXED_NUM_RR
            NUM_RR_AVERAGED=10; %this is in SECONDS!!! at least the last NUM_RR_AVERAGED seconds will be used to compute meanRR
        end
        
        % Circular buffers
        % ----------------
        BUFFER_SIZE= 5*FS;  %1250 = original value when FS = 250!!! %20000/32;%ceil((15000+16000)/32); %15000+16000 for 8000Hz, 938+1000 for 500Hz, (formula: ceil(60/HRmin*FS)+floor(60/HRmax*FS)), BUFFER_SIZE must be strictly higher than waveletSize+numberNewSamples, and higher than RRmax+numNewSamples-1 (see peakDetect.m)
        % raw data
        rawData=zeros(8,BUFFER_SIZE*(8000/FS));
        respData=zeros(1,BUFFER_SIZE*(8000/FS));
        % raw reference data
        refData=zeros(1,BUFFER_SIZE*(8000/FS));
        refrespData=zeros(1,BUFFER_SIZE*(8000/FS));
        % downsampled data and reference data
        data=zeros(8,BUFFER_SIZE);
        dataREF=zeros(1,BUFFER_SIZE);
        if RESP
            dataResp=zeros(1,BUFFER_SIZE);
            dataRespREF=zeros(1,BUFFER_SIZE);
        end
        % filtered reference data
        filteredDataREF=zeros(1,BUFFER_SIZE);
        filteredDataREF1=zeros(1,BUFFER_SIZE);
        % wavelet transformed reference data
        transformedDataREF=zeros(1,BUFFER_SIZE);
        
        % Filters
        % -------
        filterCoefficients=determineFilters(FS, ADULT); %frequency filter, temporal filter, wavelet filter,...
        filterCoefficients.iir_low1=initializeFilter(8000,1000,8,'lowpass');
        filterCoefficients.iir_low2=initializeFilter(8000,100,8,'lowpass');
        filterCoefficients.iir_low3=initializeFilter(8000,40,8,'lowpass');
        filterCoefficients.iir_high=initializeFilter(8000,100,8,'highpass');
        filterCoefficientsREF=determineFilters(FS, ADULT);
        
        % Notch filter 24Hz (@250Hz)
        % -----------------
        %         wo = 24/(FS/2); bw = wo/35;
        %         [bnotch,anotch] = iirnotch(wo,bw,-35);
        
        
        % IIR filters
        % ----------
        first=cell(1,8);
        overlap=4000;
        delay=4000;
        
        % Signal resampling
        % -----------------
        extraSamples=0;
        
        % Peak detection
        % --------------
        peakDetectParametersREF=struct('RRmin',floor(60/HR(2)*FS),'RRmax',ceil(60/HR(1)*FS),'beta',1.7,'alpha',1/3,'segmentLength',round(FS/3),'numSamplesAvailable',0,'numSamplesNeeded',FS,'threshold',[],'prevThreshold',[],'Nl',0.5,'iter',0,'previousPeak',1.2);
        peakLocationREF=[];
        O_HR_REF=0;
        
        % Averaged Heart Rate
        % -------------------
        if FIXED_NUM_RR
            RRintervalsREF=zeros(1,NUM_RR_AVERAGED);
        else
            RRintervalsREF=[];
        end
        meanRRREF=0;
        
        % Circular buffers
        % ----------------
        
        % filtered data
        filteredData=zeros(8,BUFFER_SIZE);
        filteredData1=zeros(8,BUFFER_SIZE);
        % wavelet transformed data
        transformedData=zeros(1,BUFFER_SIZE);
        
        % VCG projected data
        %projectedData=zeros(1,BUFFER_SIZE);
        % VCG projected data after Kalman filtering
        %cleanProjData=zeros(1,BUFFER_SIZE);
        
        % Channel selection
        % -----------------
        channelsToUseINIT=[1 2 3 4 5 6 7 8];
        if CHAN_SELECT==1
            proximity_matrix=[0 1 2 3 3 2 1 1;
                1 0 1 2 3 3 2 1;
                2 1 0 1 2 3 3 1 ;
                3 2 1 0 1 2 3 1 ;
                3 3 2 1 0 1 2 1 ;
                2 3 3 2 1 0 1 1 ;
                1 2 3 3 2 1 0 1 ;
                1 1 1 1 1 1 1 1 ;
                ];
        end
        
        % Reliability indicator
        % ---------------------
        old_ampli=zeros(8,3);
    end
    
    if PEAK_RESET
        %Peak detection
        %--------------
        peakDetectParameters=struct('RRmin',floor(60/HR(2)*FS),'RRmax',ceil(60/HR(1)*FS),'beta',1.7,'alpha',1/3,'segmentLength',round(FS/3),'numSamplesAvailable',0,'numSamplesNeeded',FS,'threshold',[],'prevThreshold',[],'Nl',0.5,'iter',0,'previousPeak',1.2);
        peakDetectParametersALT=struct('RRmin',floor(60/HR(2)*FS),'RRmax',ceil(60/HR(1)*FS),'beta',1.7,'alpha',1/3,'segmentLength',round(FS/3),'numSamplesAvailable',0,'numSamplesNeeded',FS,'threshold',[],'prevThreshold',[],'Nl',0.5,'iter',0,'previousPeak',1.2);
        peakLocation=[];
        if ALTERNATIVE_PEAK_DETECT
            transformedDataALT=zeros(1,BUFFER_SIZE); %AS circular buffer
            peakLocationALT=[]; %when using rel indicator OR! when detecting peaks on another signal (fusion eg)
        end
        
        % Averaged Heart Rate
        % -------------------
        if FIXED_NUM_RR
            RRintervals=zeros(1,NUM_RR_AVERAGED); %circular buffer memorize last X RR intervals
            RRintervals2=zeros(1,NUM_RR_AVERAGED);
        else
            RRintervals=[];
            RRintervals2=[];
        end
        meanRR=NaN;
        O_HR=0;
        count_holdRR=0;
    end
    
    % Kalman
    % ------
    kalmanFilterParameters=struct('N',10,'modeMeasurementNoise','auto','methodMeasurementNoise','vcg','Sigma',[],'rho',[],'psi',[],'filter',zeros(2,5),'complexStart',0.4,'complexElongation',1.4,'fadingGradientAlpha',0.2);
    [kalmanFilterParameters.filter(1,:),kalmanFilterParameters.filter(2,:)]=butter(4,2*70/FS,'low');
    kalmanFilterParameters.Sigma=zeros(8,kalmanFilterParameters.N); kalmanFilterParameters.rho=zeros(kalmanFilterParameters.N,8);
    kalmanFilterParameters.psi=zeros(8,1);
    cleanECG=struct('data',zeros(8,BUFFER_SIZE),'lastPeakLocation',[],'lastWrittenSample',[],'firstIteration','yes');
    peakToFilter=[];
    RRintervalsToFilter=[];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   OFFLINE MODE: save full signals, display graphs, create video    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if DO_INIT
        if SAVE
            % The signals to be saved (=processed windows concatenation) for later analysis
            filteredECG=[]; % fs=250Hz; matrix of 9xT where the first 8 lines are the bandpass filtered and downsampled raw cECG signals and the 9th line is the bandpass filtered and downsampled reference ECG signal
            projectedECG=[]; % fs=250Hz; the projected ECG after VCG computation (corresponds to + or - Einthoven leadII)
            cleanProjectedECG=[]; %fs = 250Hz, the projected ECG after Klaman filtering and VCG computation (corresponds to + or - Einthoven leadII)
            %vcgECG=[]; % fs=250; matrix of 2xT containing the 2D VCG components (corresponds to + or - Einthoven leadI (line1) and leadIII (line2))
            %%tempDataECG=[]; %AS
            %%tempDataECG2=[]; %AS
            referenceECG=[]; % fs=250Hz; reference ECG signal, no filtering, not synchro
            safeAmpli=[]; %fs=8000/WINDOW_SIZE Hz; amplitude of the 1kHz signal (per processing window)
            safeSelect=[]; %fs=8000/WINDOW_SIZE Hz; selected channels (per processing window)
            %safeRawData=[]; %fs=8000; 8 raw ECG signals, no filtering, no downsampling => requires a lot of memory
            %safeResp=[]; % fs=250Hz, respiration isgnal
            %safeRespREF=[]; % fs=250Hz, reference respiration signal (not synchronized!!)
            safeMotion1kHz=[]; % fs=8000/WINDOW_SIZE Hz; motion level extracted from 1kHz
            safeMotionLevel=[]; % fs=8000/WINDOW_SIZE Hz; motion level extracted from camera files (per processing window)
            safeMeanRR=[];
            safeMeanRR2=[];
            safeMeanRRREF=[];
            Rel_Indicator=[]; %fs=8000/WINDOW_SIZE Hz; reliability indicator based on heart rate and sensor coverage (per processing window)
        end
        
        if GRAPH | SAVE | VIDEO
            RRpositions=[]; %position of the QRS complexes, expressed in number of samples (fs=250Hz) from the beginning of the recording (offset included!!), take diff(RRpositions) to get the value of the R-R intervals
            RRpositionsALT=[]; %AS
            RRpositionsREF=[];
        end
        
        if GRAPH
            % Dynamic graph
            time=[1:BUFFER_SIZE]./FS;
            fig=figure;
            set(fig,'doublebuffer','on')
            ax(1)=subplot(4,1,1); plot(time, dataREF, '.r'); title('Reference ECG')
            ax(2)=subplot(4,1,2); plot(time, data(1,:), '.'); title('LeadII projected data')
            ax(3)=subplot(4,1,3); plot(time, data(1,:), '.'); title('Wavelet transform and detected peaks')
            ax(4)=subplot(4,1,4); plot(time, cleanECG.data(1,:),'.'); title('Cleaned leadII projected data')
            xlabel('time (seconds)')
            [hr, el, abno, rel] = template_dynamic_graph(fig);
        end
        
        if VIDEO
            % Video initialization
            if RESP
                [figVid, ha, el] = template_video_resp;
            else
                [figVid, ha, el] = template_video;
            end
            time=[1:BUFFER_SIZE]./FS;
            plot(ha(1),time, dataREF, '.r');
            vidObj = VideoWriter('Video_neonat.avi');
            vidObj.FrameRate=8000/numberNewSamples;
            vidObj.Quality=80; %quality: 90
            open(vidObj);
        end
    end
    % Reset
    % -----
    countBeforeReset=0;
    DO_INIT=0;
    PEAK_RESET=0;
    KALMAN_RESET=0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         INFINITE LOOP                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
%-----------------------Process inputs------------------------------------%
%-------------------------------------------------------------------------%

% Determine the VCG projection angle knowing the body position
% Lying supine:0 , prone:1  , Left side:2,  Right side:3 , Automated:4
if Manual==1
    updatePreferentialDirection=0;
    if Body==1 %chest
        preferentialDirection=-123/180*pi;
    elseif Body==0 %back
        preferentialDirection=-57/180*pi;
    else
        updatePreferentialDirection=1;
    end
else
    % Allow automated body position detection
    updatePreferentialDirection=1;
end
% Output
O_Body=Body;
% Get new data
newrefData=ECG1;
newData=[SEN1;SEN2;SEN3;SEN4;SEN5;SEN6;SEN7;SEN8];

numberNewSamples=size(newData,2); %the window size can vary over time (especially for online demo
if ONLINE_DEMO & RESP
    referenceRespi=ECG1; %the reference respiration is not yet available for online demo
end
%Update circular buffer with new data samples
rawData=[rawData(:,numberNewSamples+1:end),newData];
refData=[refData(1,numberNewSamples+1:end),newrefData];
if RESP
    refrespData=[refrespData(1,numberNewSamples+1:end),referenceRespi];
end

%-------------------------------------------------------------------------%
%------------Perform channel selection using 1kHz (SENSOR COVERAGE)-------%
%-------------------------------------------------------------------------%

% PERFORM HERE THE CHANNEL SELECTION ------->>>>>>>>
% ----------------------------------
% Use whatever strategy you want using as input the rawData matrix and as output
% the list of channels to use (e.g. channelsToUse=[1 2 4 6 8];)

if CHAN_SELECT==0 || CHAN_SELECT==3
    % Aline's strategy:
    if ~exist('indexMax','var')
        [maxi,indexMax]=max(rawData(:,end-160+1:end),[],2);
        indexMax=indexMax+BUFFER_SIZE*(8000/FS)-160;
    else
        indexMax=indexMax-mod(numberNewSamples,160);
        if indexMax<160 | rawData(indexMax)==0
            [maxi,indexMax]=max(rawData(:,end-160+1:end),[],2);
            indexMax=indexMax+BUFFER_SIZE*(8000/FS)-160;
        end
    end
    [maxi,i1]=max(rawData(:,indexMax-3:min(indexMax+3,BUFFER_SIZE*(8000/FS))),[],2);
    [mini,i2]=min(rawData(:,indexMax:min(indexMax+5,BUFFER_SIZE*(8000/FS))),[],2);
    ampli=maxi-mini; % needed for the rel indicator!
    channelsToUse=channelsToUseINIT(ampli<=1.25*median(ampli)&ampli<=THRESHOLD&ampli>MIN_THRESHOLD);%1.5*median before % not needed for CHAN_SELECT==3!
    
elseif CHAN_SELECT==1
    % Louis's strategy:
    % Simply high pass filtering here
    fs=8000;
    fc_hpf=800;
    %Sort out if there are channels less than a certain value
    %useful=find(max(rawData,[],2)>1000);
    useful=channelsToUseINIT(max(rawData,[],2)>1000); %AS
    if ~isempty(useful)
        fNorm = fc_hpf / (fs/2);
        [d,c] = butter(8, fNorm, 'high');
        rawData1= (filtfilt(d, c, (rawData(useful,:))'))';
        chunker=12;
        [loc,z]=outspikeloc(rawData1(:,100:end));
        mranvec=meanerdetail(z', chunker);
        [rt,chorder]=sort(mranvec);
        bestchannel=chorder(1);
        cors=corrcoef(rawData(:,100:5:end)');
        [ty,yu]=sort(cors(:,bestchannel));
        %Here do something about proximity:
        outs= flip_proximity(yu,proximity_matrix);
        channelsToUse=outs(1:3);
        %  channelsToUse=flipud(yu(end-2:end));
        % channelsToUse=chorder(end-2:end)
        %         if SAVE
        %             channeltrail=[channeltrail;channelsToUse];
        %         end
        [maxi,i1]=max(rawData,[],2);
        [mini,i2]=min(rawData,[],2);
        %         %
        %         % %      figure;plot(rawData(:,indexMax-3:min(indexMax+3,BUFFER_SIZE*(8000/FS))))
        %         % %      legend('1','2','3','4','5','6','7','8'
        ampli=maxi-mini;
    else
        channelsToUse=[2 4 8];
    end
    % elseif CHAN_SELECT==2
    % Mohammed's strategy:
    % take always 2 best channels
    % channelsToUse=[...];
    % else
    %     channelsToUse=channelsToUseINIT;
end

% Mohammed/Louis: discard here any broken sensor if the MIN_THRESHOLD is not enough:
% e.g. 5 is broken:
% channelsToUse=channelsToUse(channelsToUse~=5);

% <<<<<<<<<----------


% PERFORM A PARALLEL SELECTION OF !2 CHANNELS MAX! (for comparison purpose)
% ------------------------------------------------
% Here you can select ONLY 2 channels, a parallel peak detection will be
% performed on channelsToUseALT(1)-channelsToUseALT(2) and the results
% written in the text file
if CHAN_SELECT~=3
    % Mohammed's strategy:
    if ALTERNATIVE_PEAK_DETECT
        if CHAN_2SELECT==0
            % Always take chan2-chan4 for peak detection
            channelsToUseALT(1)=2;
            channelsToUseALT(2)=4;
        elseif CHAN_2SELECT==1
            % Always take the 2 best coupled channels for peak detection
            if CHAN_SELECT==0
                % Aline's selection
                [junk,channelsToUseALT]=sort(ampli);
            elseif CHAN_SELECT==1
                % Louis'selection
                channelsToUseALT(1)=channelsToUse(1);
                channelsToUseALT(2)=channelsToUse(2);
            else
                channelsToUseALT=channelsToUse;
            end
        elseif CHAN_2SELECT==2
            % Always take the 2 best non-adjacent coupled channels (otherwise,chan 2 - chan 4)
            if (channelsToUse(1)-channelsToUse(2))<=1
                if (channelsToUse(1)-channelsToUse(3))<=1
                    if (channelsToUse(2)-channelsToUse(3))<=1
                        channelsToUseALT=[2 4];
                    else
                        channelsToUseALT=[channelsToUse(2) channelsToUse(3)];
                    end
                else
                    channelsToUseALT=[channelsToUse(1) channelsToUse(3)];
                end
            else
                channelsToUseALT=channelsToUse(1:2);
            end
        end
    end
    
    
    % FORCE THE SELECTION OF AT LEAST 3 CHANNELS:
    % -------------------------------------------
    % This makes sure that
    if  MIN3CHANNELS & numel(channelsToUse)<=2
        %disp(['Step1before: ' num2str(channelsToUse)])
        channelsToUse=channelsToUseINIT(ampli<=1.5*median(ampli)&ampli<=THRESHOLD&ampli>MIN_THRESHOLD);%1300000);
        %disp(['Step1after: ' num2str(channelsToUse)])
        if numel(channelsToUse)<=2
            channelsToUse=channelsToUseINIT(ampli<=1.5*median(ampli)&ampli<=1.5*THRESHOLD&ampli>MIN_THRESHOLD);
            %disp(['Step2after: ' num2str(channelsToUse)])
            %TO DO: here you should actually take the diff between the two channels
            %and NOT compute the VCG
            if numel(channelsToUse)<=2
                channelsToUse=channelsToUseINIT(ampli<=1.5*median(ampli)&ampli<=2*THRESHOLD&ampli>MIN_THRESHOLD);
                %disp(['Step3after: ' num2str(channelsToUse)])
                if numel(channelsToUse)<2
                    if VERBOSE
                        disp('0 or 1 sensor selected: impossible to find 3 good coupled channels')
                    end
                elseif numel(channelsToUse)==2
                    if VERBOSE
                        disp(['2 SENSORS SELECTED: 2 best coupled channels were selected because impossible to find 3 good coupled channels:' num2str(channelsToUse)])
                    end
                end
            end
        end
    end
    
end



%-------------------------------------------------------------------------%
%-------------------Respiration extraction (Louis)------------------------%
%-------------------------------------------------------------------------%

if RESP
    % AS: This needs to be optimized and verified
    k=8;
    %for k=1:8
    raw{k}=rawData(k,end-numberNewSamples-overlap-delay+1:end);
    env1{k}=hilbertdemod(raw{k},filterCoefficients,first{k},overlap+delay);
    % to avoid jumps in case of channel change
    env1{k}=env1{k}-mean(env1{k},1);
    %
    if isempty(first{k})
        first{k}=0;
        [prevs1{k},s1{k}]=iirFilter([],env1{k},filterCoefficients.iir_low3,'first',overlap+delay);
    else
        [prevs1{k},s1{k}]=iirFilter(prevs1{k},env1{k},filterCoefficients.iir_low3,'',overlap+delay);
    end
    s1{k}=smooth(s1{k},2000);
    %end
    
    %     %Now for channel selection, we're choosing the best one here:
    %     channelssorted=channelrangecorr(raw,s1,1,length(s1{1}),2);
    %     %This step contains the resp signal for this window
    %     windowedrresult=windowedresp{channelssorted(1)}(overlap+1:end-delay); %should remove the 600 last (delay compensation?!)?!
    
    windowedrresult=s1{k}(overlap+1:end-delay); %force channel 4
    
    respData=[respData(1,numberNewSamples+1:end),windowedrresult'];
    respData=smooth(respData,1000)';
    
    
    
    %-------------------------------------------------------------------------%
    %-------------------Filter and downsample data ---------------------------%
    %-------------------------------------------------------------------------%
    
    % if RESP:
    [dataOut,extraSamples,numberNewSamples]=signalResampling([rawData;refData;respData;refrespData],[data;dataREF;dataResp;dataRespREF],filterCoefficients,numberNewSamples,extraSamples,(8000/FS)); %32=8000/250
    data=dataOut(1:8,:);
    dataREF=dataOut(9,:);
    dataResp=dataOut(10,:);
    dataRespREF=dataOut(11,:);
    if SAVE
        %safeResp=[safeResp, dataResp];
        %safeRespREF=[safeRespREF, dataRespREF];
        safeMotion1kHz=[safeMotion1kHz, var(dataResp(end-numberNewSamples+1:end))]; %mean(dataResp(end-numberNewSamples+1:end),2)
    end
else
    [dataOut,extraSamples,numberNewSamples]=signalResampling([rawData;refData],[data;dataREF],filterCoefficients,numberNewSamples,extraSamples,(8000/FS)); %32=8000/250
    data=dataOut(1:8,:);
    dataREF=dataOut(9,:);
end

if SAVE
    referenceECG=[referenceECG,dataREF(end-numberNewSamples+1:end)]; %NOT filtered so NOT synchronized
end

%---------------------------------------------------------------------%
%-------------------Apply temporal filter (TBD)-----------------------%
%---------------------------------------------------------------------%

%Apply temporal filter on data to compensate for the transfer function
%of the capacitive sensors (although in analysis the temporal filter
%was applied/derived from a spatially combined ECG, it is a linear
%filter and can hence also be applied in this stage already)
%filteredData=temporalFilter(data,filteredData,filterCoefficients.sensorTransferFunction,numberNewSamples);
%Note: data WILL BE DELAYED BY 3/250=12ms if FS = 250Hz

%---------------------------------------------------------------------%
%------------------------Apply band-pass filter-----------------------%
%---------------------------------------------------------------------%

filteredData1=bandpassFiltering(filteredData1, data, filterCoefficients, numberNewSamples,1);
filteredDataREF1=bandpassFiltering(filteredDataREF1, dataREF, filterCoefficients, numberNewSamples,1);

%---------------------------------------------------------------------%
%------------------------Apply NOTCH filter (24Hz) -------------------%
%---------------------------------------------------------------------%

if NOTCH24
    filteredData=notchFiltering(filteredData, filteredData1, filterCoefficients, numberNewSamples);
    filteredDataREF=notchFiltering(filteredDataREF, filteredDataREF1, filterCoefficients, numberNewSamples);
    % Perform ref-based channel selection
else
    filteredData=filteredData1;
    filteredDataREF=filteredDataREF1;
end

%---------------------------------------------------------------------%
%-----------------------Channel selection based on REF----------------%
%---------------------------------------------------------------------%

if CHAN_SELECT==3
    [S,ord]= spotdiff(filteredData,filteredDataREF);
    channelsToUse=ord;
    if ALTERNATIVE_PEAK_DETECT==1
        channelsToUseALT=ord;
    end
end

%-------------------------------------------------------------------------%
%--------------Save some variables for later offline analysis ------------%
%-------------------------------------------------------------------------%

if ONLINE_DEMO | SAVE | VIDEO | GRAPH
    % This is used to display the flower colors...
    sensor_coverage=3*ones(1,8);
    if ADULT
        sensor_coverage(channelsToUse)=1; %1=green, 2=orange, 3=red
    else %the neonatal nattress has a different senor configuration (turned 180degrees)
        shiftedChan=channelsToUse+3;
        shiftedChan(shiftedChan==8)=1;
        shiftedChan(shiftedChan==11)=8;
        shiftedChan(shiftedChan>8)=shiftedChan(shiftedChan>8)-7;
        sensor_coverage(shiftedChan)=1;
    end
end


if SAVE
    safeAmpli=[safeAmpli, ampli];
    %safeRawData=[safeRawData, newData]; %newrefData
    safeSelect=[safeSelect, sensor_coverage'];
end


%-------------------------------------------------------------------------%
%-----------Decide to do processing and display results of not------------%
%-------------------------------------------------------------------------%
% If not enough channels are reliable, the processing is paused and zeros
% are displayed. If this bad condition remains for more than X counts
% (if X=3, 3*0.5 = 1.5 seconds) then the algorithm is restarted
if numel(channelsToUse)<=1
    % Process the reference only
    transformedDataREF=signalTransformation2(filteredDataREF,transformedDataREF,filterCoefficientsREF,numberNewSamples);
    [peakDetectParametersREF,peakLocationREF,newRRintervalsREF]=peakDetect2(transformedDataREF,peakDetectParametersREF,numberNewSamples,FS);
    if numel(newRRintervalsREF)>=1
        if FIXED_NUM_RR
            %OPTION1: compute the mean on the last "NUM_RR_AVERAGED" RR intervals AS
            RRintervalsREF=[RRintervalsREF(numel(newRRintervalsREF)+1:end), newRRintervalsREF];
        else
            %OPTION2: compute the mean on the last 10 seconds AS
            RRintervalsREF=[RRintervalsREF, newRRintervalsREF];
            while sum(RRintervalsREF(2:end))>=NUM_RR_AVERAGED*FS %10seconds
                RRintervalsREF(1)=[]; %remove oldest RRinterval
            end
        end
        %Average HR over numel(RRintervals) last values
        meanRRREF=mean(RRintervalsREF(RRintervalsREF~=0),2); %AS: TODO: progressive update of the mean
        O_HR_REF=round(60*FS/meanRRREF);
    end
    % Outputs
    Abnormality='Lost contact';
    O_ref_ECG=filteredDataREF(1,end-numberNewSamples+1:end);
    proc_ECG_lead2=zeros(1,numberNewSamples);
    proc_ECG_lead1=zeros(1,numberNewSamples);
    proc_ECG_lead3=zeros(1,numberNewSamples);
    proc_ECG_lead2_blue=zeros(1,numberNewSamples);
    proc_ECG_lead1_blue=zeros(1,numberNewSamples);
    proc_ECG_lead3_blue=zeros(1,numberNewSamples);
    if numel(RRintervals)>=3
        meanRR=mean(RRintervals(RRintervals~=0),2); %AS: TODO: progressive update of the mean &do it over 10seconds (while cumsum < 10*250, add it)
    end
    O_HR2=0;
    Rel_Indc=1;
    cleanProjData=NaN+zeros(1,BUFFER_SIZE);
    projectedData=NaN+zeros(1,BUFFER_SIZE);
    % Reset
    countBeforeReset=countBeforeReset+1;
    if countBeforeReset>=NUMBER_BEFORE_RESET
        Abnormality='RESET';
        KALMAN_RESET=1; PEAK_RESET=1;
    end
    % Reliability indicator
    previousIsBad=1;
    
    % mean HR
    count_holdRR=count_holdRR+1;
    if count_holdRR*WINDOW_SIZE/8000 > ceil(60/HR(1))  %holding value for more than the longest RR interval possible ---> not possible ! --> Reset!
        O_HR=NaN;
        if SAVE
            safeMeanRR=[safeMeanRR, NaN];
        end
        KALMAN_RESET=1; PEAK_RESET=1;
    else
        O_HR=round(60*FS/meanRR); %hold
        if SAVE
            safeMeanRR=[safeMeanRR, meanRR]; %AS:safeMeanRR(end)
        end
    end
    if SAVE
        safeMeanRRREF=[safeMeanRRREF, meanRRREF]; %safeMeanRRREF(end)];
        safeMeanRR2=[safeMeanRR2, meanRR2];
    end
else
    Abnormality='none';
    countBeforeReset=0;
    
    
    %---------------------------------------------------------------------%
    %-------------------Take bipolar leads--------------------------------%
    %---------------------------------------------------------------------%
    
    %Use mean of all signals as reference. With approach of combining signals into VCG, the electrode in the middle(#8) will hardly be used as it is so close to the reference (the mean).
    %Sanity check required to see whether any of the signals can be discarded
    %and excluded from the determination of the mean
    %channelsToUse=sanityCheckForMean(data); %Decision based on the last 2 seconds of data
    filteredData(:,BUFFER_SIZE-numberNewSamples+1:BUFFER_SIZE)=filteredData(:,BUFFER_SIZE-numberNewSamples+1:BUFFER_SIZE)-repmat(mean(filteredData(channelsToUse,BUFFER_SIZE-numberNewSamples+1:BUFFER_SIZE)),size(filteredData,1),1);
    %---------------------------------------------------------------------%
    %-------------------Find transformation matrixes----------------------%
    %---------------------------------------------------------------------%
    
    if numel(channelsToUse) > 2 || ONLINE_DEMO || KALMAN
        if ~ADULT
            [vcg2ecgMatrix,transformationMatrix]=findTransformationMatrix('neonatalV2',channelsToUse);
        else
            [vcg2ecgMatrix,transformationMatrix]=findTransformationMatrix('adult',channelsToUse); %TBD: update adult array distances
        end
        
        %---------------------------------------------------------------------%
        %-------------------Compute VCG (TBD)---------------------------------%
        %---------------------------------------------------------------------%
        
        % Find preferential direction for VCG projection that enhance the
        % R-peaks amplitude
        if updatePreferentialDirection==1
            %Sanity check on all incoming data channels
            %[channelsToUse,vcg]=sanityCheck(filteredData,vcg2ecgMatrix,transformationMatrix,channelsToUse);
            %%AS: to be studied
            %TBD: Automatically find preferential direction in VCG that yields maximum QRS amplitude (based on
            %the last BUFFER_SIZE/FS seconds of data)
            %The preferential direction for the data will be determined based on the
            %last 5 seconds of data. Newly arriving data will be projected on the last
            %known direction until enough new data is there to compute a new direction.
            %[preferentialDirection,oldPreferentialDirections]=findPreferentialDirection(vcg,oldPreferentialDirections);
            preferentialDirection=57/180*pi;
        end
        vcg=transformationMatrix*filteredData(channelsToUse,:);  %AS: TO DO: only on new samples! + use sanity check for improved channel selection
    end
    %---------------------------------------------------------------------%
    %-------------------Project data on 3 Einthoven leads-----------------%
    %---------------------------------------------------------------------%
    if numel(channelsToUse) > 2
        projectedData=[cos(preferentialDirection) sin(preferentialDirection)]*vcg; %will be Inf and -Inf if channelsToUse==2
    else %the tranformationMatrix is NaN if only 2 channels ---> no vcg, no reprojection
        projectedData=filteredData(channelsToUse(1),:)-filteredData(channelsToUse(2),:);
    end
    
    % Outputs of the online demo (first and second graphs)
    % --------------------------
    if ONLINE_DEMO
        % Reference ecg
        O_ref_ECG=filteredDataREF(1,end-numberNewSamples+1:end);
        % LEADII
        proc_ECG_lead2=-projectedData(1,end-numberNewSamples+1:end); %inverted input signals (-)
        % LEADI
        if Body==0 %back
            proc_ECG_lead1=-vcg(1,end-numberNewSamples+1:end); %inverted input signals (-)
        else %chest etc
            proc_ECG_lead1=vcg(1,end-numberNewSamples+1:end);  %inverted input signals (-) and lead1 = opposite X axis (-)
        end
        % LEADIII
        proc_ECG_lead3=vcg(2,end-numberNewSamples+1:end); %inverted input signals (-) and lead3 = bottom to top =>different from Y axis (-)
    end
    
    if SAVE
        %vcgECG=[vcgECG,vcg(:,end-numberNewSamples+1:end)];
    end
    
    %---------------------------------------------------------------------%
    %-------------------Wavelet transform---------------------------------%
    %---------------------------------------------------------------------%
    
    %Transform data (using wavelet convolution) to facilitate peak
    %detection. Note that the last waveletSize/2 samples are not
    %defined (corresponds to a delay in the real-time output)
    transformedDataREF=signalTransformation2(filteredDataREF,transformedDataREF,filterCoefficients,numberNewSamples);
    transformedData=signalTransformation2(projectedData,transformedData,filterCoefficients,numberNewSamples);
    
    if ALTERNATIVE_PEAK_DETECT
        %transformedDataALT=signalTransformation2(filteredData(2,:)-filteredData(4,:),transformedDataALT,filterCoefficients,numberNewSamples);%AS
        transformedDataALT=signalTransformation2(filteredData(channelsToUseALT(1),:)-filteredData(channelsToUseALT(2),:),transformedDataALT,filterCoefficients,numberNewSamples);%AS
    end
    
    %---------------------------------------------------------------------%
    %-------------------Peak detection------------------------------------%
    %---------------------------------------------------------------------%
    %Peak detection itself is used to find the heartrate and enable
    %further enhancement of the ECG (see HR-based Kalman filter)
    %NB:AS: <peakLocation> is given in samples from start of the vector <transformedData>
    %<peakLocation> is empty or is a line vector of a few peaks
    %NB:AS: problem with iter2 and the update of prevThreshold
    
    [peakDetectParametersREF,peakLocationREF,newRRintervalsREF]=peakDetect2(transformedDataREF,peakDetectParametersREF,numberNewSamples,FS);
    [peakDetectParameters,peakLocation,newRRintervals]=peakDetect2(transformedData,peakDetectParameters,numberNewSamples,FS);
    if ALTERNATIVE_PEAK_DETECT
        [peakDetectParametersALT,peakLocationALT,junk]=peakDetect2(transformedDataALT,peakDetectParametersALT,numberNewSamples,FS); %AS
    end
    %AS: perform peak fusion or alternative peak detection
    
    %---------------------------------------------------------------------%
    %-------------------Averaged heart rate-------------------------------%
    %---------------------------------------------------------------------%
    %Update buffer with RR intervals (in number of samples)
    if numel(newRRintervals)>=1
        count_holdRR=0;
        %RR intervals from capacitive sensors
        if FIXED_NUM_RR
            %OPTION1: compute the mean on the last "NUM_RR_AVERAGED" RR intervals AS
            RRintervals=[RRintervals(numel(newRRintervals)+1:end), newRRintervals];
        else
            %OPTION2: compute the mean on the last 10 seconds AS
            RRintervals=[RRintervals, newRRintervals];
            while sum(RRintervals(2:end))>=NUM_RR_AVERAGED*FS %10seconds
                RRintervals(1)=[]; %remove oldest RRinterval
            end
        end
        %Average HR over numel(RRintervals) last values
        meanRR=mean(RRintervals(RRintervals~=0),2); %TODO: progressive update of the mean
        O_HR=round(60*FS/meanRR);
        if SAVE
            safeMeanRR=[safeMeanRR, meanRR];
        end
    else %hold or output NaN
        count_holdRR=count_holdRR+1;
        if count_holdRR*WINDOW_SIZE/8000 > ceil(60/HR(1))
            O_HR=NaN;
            if SAVE
                safeMeanRR=[safeMeanRR, NaN];
            end
            KALMAN_RESET=1; PEAK_RESET=1;
        else
            O_HR=round(60*FS/meanRR); %hold
            if SAVE
                safeMeanRR=[safeMeanRR, meanRR]; %AS:safeMeanRR(end)
            end
            
        end
    end
    if numel(newRRintervalsREF)>=1
        %RR intervals from reference signal
        if FIXED_NUM_RR
            %OPTION1: compute the mean on the last "NUM_RR_AVERAGED" RR intervals AS
            RRintervalsREF=[RRintervalsREF(numel(newRRintervalsREF)+1:end), newRRintervalsREF];
        else
            %OPTION2: compute the mean on the last 10 seconds AS
            RRintervalsREF=[RRintervalsREF, newRRintervalsREF];
            while sum(RRintervalsREF(2:end))>=NUM_RR_AVERAGED*FS %10seconds
                RRintervalsREF(1)=[]; %remove oldest RRinterval
            end
        end
        %Average HR over numel(RRintervals) last values
        meanRRREF=mean(RRintervalsREF(RRintervalsREF~=0),2); %TODO: progressive update of the mean
        O_HR_REF=round(60*FS/meanRRREF);
        if SAVE
            safeMeanRRREF=[safeMeanRRREF, meanRRREF];
        end
    else
        if SAVE
            if isempty(safeMeanRRREF)
                safeMeanRRREF=0;
            else
                safeMeanRRREF=[safeMeanRRREF, meanRRREF]; %safeMeanRRREF(end)];
            end
        end
    end
    
    
    %---------------------------------------------------------------------%
    %------------------------Reliability indicator------------------------%
    %---------------------------------------------------------------------%
    
    % Quick and dirty implementation of a reliability indicator based on
    % heart rate and sensor coverage
    % Goal: call the Kalman filter only when signal is reliable!
    old_ampli=[old_ampli(:,2:end) ampli];
    if numel(channelsToUse)>=2 & sum((old_ampli(channelsToUse,end)-old_ampli(channelsToUse,end-1))<=old_ampli(channelsToUse,end-1)/5)>=min(3,numel(channelsToUse)-1)
        if isempty(newRRintervals)
            cumu=cumu+peakDetectParameters.RRmin+abs(peakDetectParameters.previousPeak-1)*FS+numberNewSamples-peakDetectParameters.numSamplesAvailable; %compute an RR interval from the last peak ntil the end of the new samples processed (where nothing was found)
            if previousIsBad~=0 | cumu> 1.1*meanRR2
                Rel_Indc=1;
                O_HR2=0;
            end
        elseif(60*FS/HR(2)<=newRRintervals & newRRintervals<=60*FS/HR(1) & newRRintervals<=1.7*meanRR2 & newRRintervals>=meanRR2/1.7)
            cumu=0;
            if FIXED_NUM_RR
                %OPTION1: compute the mean on the last "NUM_RR_AVERAGED" RR intervals AS
                RRintervals2=[RRintervals2(numel(newRRintervals)+1:end), newRRintervals];
            else
                %OPTION2: compute the mean on 10 seconds (not necessarly the last)AS
                RRintervals2=[RRintervals2, newRRintervals];
                while sum(RRintervals2(2:end))>=NUM_RR_AVERAGED*FS %10seconds
                    RRintervals2(1)=[]; %remove oldest RRintervals
                end
            end
            if numel(RRintervals2(RRintervals2~=0))>=3
                %                 if previousIsBad==1
                %                     % one good RR interval is needed to reactivate the trust
                %                     previousIsBad=0.5; %0.5 %indicates that if the next RR in the bounds, it will reactivate trust
                %                     Rel_Indc=1;
                %                     O_HR2=0;
                %                 else
                %the RR intervals reactivates trust
                previousIsBad=0;
                %peakLocation2=peakLocation; %AS change me peakLocationALT
                Rel_Indc=0;
                ind=find(RRintervals2~=0,1);
                if ~isempty(ind)&ind~=1
                    meanRR2=mean(RRintervals2(ind+2:end),2); %TODO: progressive update of the mean
                    O_HR2=round(60*FS/meanRR2);
                else
                    meanRR2=mean(RRintervals2,2); %TODO: progressive update of the mean
                    O_HR2=round(60*FS/meanRR2);
                end
                %                 end
            else %not enough past RRintervals
                Rel_Indc=1; %bad
                previousIsBad=1;
                O_HR2=0;
                Abnormality='RR';
            end
        else % out of bound RR intervals
            Rel_Indc=1; %bad
            previousIsBad=1;
            O_HR2=0;
            Abnormality='RR';
        end
    else % bad coupling %the bad coupling is actually a prediction about 0.8 s later (prediction made before the filtering bank => not yet delayed)
        cumu=0;
        Rel_Indc=1;
        previousIsBad=1;
        O_HR2=0;
        Abnormality='Bad coupling';
    end
    
    if SAVE
        if O_HR2==0
            safeMeanRR2=[safeMeanRR2, NaN];
        else
            safeMeanRR2=[safeMeanRR2, meanRR2];
        end
    end
    %---------------------------------------------------------------------%
    %-------------------Signal enhancement (Kalman filter)----------------%
    %---------------------------------------------------------------------%
    
    %AS: implementation could be optimized and simplified
    %if Rel_Indc==0 %AS
    if KALMAN
        %Filter ECG based on quasi-periodicity
        cleanECG.data=[cleanECG.data(:,numberNewSamples+1:end) zeros(8,numberNewSamples)];
        cleanECG.lastPeakLocation=cleanECG.lastPeakLocation-numberNewSamples;
        cleanECG.lastWrittenSample=cleanECG.lastWrittenSample-numberNewSamples;
        if isempty(peakToFilter)==0
            peakToFilter=peakToFilter-numberNewSamples;
        end
        peakToFilter=[peakToFilter peakLocation];
        RRintervalsToFilter=[RRintervalsToFilter newRRintervals];
        
        if isempty(RRintervalsToFilter)==0
            h=1;
            while h<=length(peakToFilter)
                %NB: use [projectedData;vcg] instead of filtered data AS:TODO
                [cleanECG,kalmanFilterParameters,status,KALMAN_RESET]=kalmanFilter(cleanECG,kalmanFilterParameters,filteredData(channelsToUse,:),peakToFilter(h),RRintervalsToFilter(h),transformationMatrix,vcg2ecgMatrix,channelsToUse,KALMAN_RESET);
                if strcmp(status,'notFiltered')==0
                    %Omit current <peakToFilter> and <instantaneousRRintervals>
                    peakToFilter=peakToFilter(h+1:end);
                    RRintervalsToFilter=RRintervalsToFilter(h+1:end);
                    %don't increment h. The discarding in peakToFilter and
                    %RRintervals has the same effect.
                else
                    h=h+1; %go to see whether, for some strange reason, the current complex can't be filtered, but the next one can. Probably this will happen without result, but this implementation is more easy to make than other implementations that check all complexes
                end
            end
        end
    end
    
    %---------------------------------------------------------------------%
    %-----------------------Project on Einthoven leads--------------------%
    %---------------------------------------------------------------------%
    
    if numel(channelsToUse) > 2
        vcg=transformationMatrix*cleanECG.data(channelsToUse,:);
        cleanProjData=[cos(preferentialDirection) sin(preferentialDirection)]*vcg;
    else %the tranformationMatrix is NaN if only 2 channels ---> no vcg, no reprojection
        cleanProjData=cleanECG.data(channelsToUse(1),:)-cleanECG.data(channelsToUse(2),:);
    end
    
    
    % Outputs of the online demo (third graph)
    % ----------------------------------------
    % Signals after Klaman smoothing
    
    if ONLINE_DEMO
        % LEADII
        proc_ECG_lead2_blue=-cleanProjData(1,end-(FS+FS/2)+1:end-FS); %inverted input signals (-)
        % LEADI
        if RESP
            proc_ECG_lead1_blue=dataResp(1,end-numberNewSamples+1:end); %output the resp signal instead of leadI
        else
            if Body==0 %back
                proc_ECG_lead1_blue=-vcg(1,end-(FS+FS/2)+1:end-FS); %inverted input signals (-)
            else %chest etc
                proc_ECG_lead1_blue=vcg(1,end-(FS+FS/2)+1:end-FS);  %inverted input signals (-) and lead1 = opposite X axis (-)
            end
        end
        % LEADIII
        %proc_ECG_lead3_blue=vcg(2,end-375+1:end-250); %inverted input signals (-) and lead3 = bottom to top =>different from Y axis (-)
        proc_ECG_lead3_blue=transformedData(1,end-numberNewSamples+1:end); %output the wavelet transform (=peak position) instead of leadIII
    end
end

%-------------------------------------------------------------------------%
%------------------------------Save signals-------------------------------%
%-------------------------------------------------------------------------%

if SAVE
    Rel_Indicator=[Rel_Indicator, Rel_Indc];
    %Rel_Indicator=[Rel_Indicator, [Rel_Indc; numberNewSamples] ]; %for variable window lengths
    toBeAdded= [filteredData(:,end-numberNewSamples+1:end); filteredDataREF(end-numberNewSamples+1:end)];
    filteredECG=[filteredECG, toBeAdded]; %the 8 sensors plus the filtered reference ECG
    projectedECG=[projectedECG,projectedData(end-numberNewSamples+1:end)];
    cleanProjectedECG=[cleanProjectedECG,cleanProjData(1,end-375-numberNewSamples+1:end-375)];
end
if SAVE | GRAPH | VIDEO
    %Save all RRpositions of the dataset (OFFLINE_MODE only)
    for i=1:numel(peakLocation)
        RRpositions=[RRpositions, peakLocation(i)+floor(sampleCounter/(8000/FS))-BUFFER_SIZE]; %take diff(RRpositions) to get the RRintervals
    end
    if ALTERNATIVE_PEAK_DETECT
        for i=1:numel(peakLocationALT) %AS change me
            RRpositionsALT=[RRpositionsALT, peakLocationALT(i)+floor(sampleCounter/(8000/FS))-BUFFER_SIZE]; %take diff(RRpositions) to get the RRintervals (-1+NNS)
        end
    end
    for i=1:numel(peakLocationREF)
        RRpositionsREF=[RRpositionsREF, peakLocationREF(i)+sampleCounter/(8000/FS)-BUFFER_SIZE]; %take diff(RRpositions) (+NNS*FSRAT-1) to get the RRintervals
    end
end

%-------------------------------------------------------------------------%
%---------------------------Update dynamic graph--------------------------%
%-------------------------------------------------------------------------%

if GRAPH
    plot(ax(1),time,filteredDataREF, 'r'); title(ax(1),'Reference ECG')
    plot(ax(2),time,-projectedData); title(ax(2),'LeadII projected data')
    plot(ax(3),time,-cleanProjData); title(ax(3),'Cleaned leadII data')
    if RESP
        plot(ax(3),time,dataRespREF); title(ax(3),'Reference Respiration')
        plot(ax(4),time,dataResp,'lineWidth',1.5,'color',[0 0.75 0.75]); title(ax(4),'Respiration');
    else
        if numel(RRpositions)>=2
            irr=diff(RRpositions);
            stairs(ax(4), RRpositions(max(2,find(RRpositions>=RRpositions(end-1)-15000)):end)./FS,60*FS./irr(max(2,find(RRpositions>=RRpositions(end-1)-15000))-1:end),'b'); title(ax(4),'Instantaneous HR'); ylim(ax(4),[45 85]);
        end
        if numel(RRpositionsALT)>=2
            axes(ax(4));
            hold on
            irr=diff(RRpositionsALT);
            stairs(ax(4), RRpositionsALT(max(2,find(RRpositionsALT>=RRpositionsALT(end-1)-15000)):end)./FS,60*FS./irr(max(2,find(RRpositionsALT>=RRpositionsALT(end-1)-15000))-1:end),'g'); title(ax(4),'Instantaneous HR'); ylim(ax(4),[45 85]);
        end
        if numel(RRpositionsREF)>=2
            axes(ax(4));
            hold on
            irr=diff(RRpositionsREF);
            stairs(ax(4), RRpositionsREF(max(2,find(RRpositionsREF>=RRpositionsREF(end-1)-15000)):end)./FS,60*FS./irr(max(2,find(RRpositionsREF>=RRpositionsREF(end-1)-15000))-1:end),'r'); title(ax(4),'Instantaneous HR'); ylim(ax(4),[45 85]);
            hold off;
        end
    end
    set(hr(1),'String',{num2str(O_HR_REF)});
    set(hr(2),'String',{num2str(O_HR)});
    set(hr(3),'String',{num2str(O_HR2)});
    %Sensor coverage
    set(el(1:8),'FaceColor',[0 0 0]);
    set(el(sensor_coverage==1),'FaceColor',[0 0.7 0.1]);
    set(el(sensor_coverage==2),'FaceColor',[0.5 0.7 0]);
    set(el(sensor_coverage==3),'FaceColor',[1 0 0]);
    
    set(abno(1),'String',Abnormality);
    if Rel_Indc==0
        set(rel(1),'String','Good','Color',[0 1 0]);
    else
        set(rel(1),'String','Bad','Color',[1 0 0]);
    end
    
    for i=1:numel(peakLocation)
        axes(ax(2));
        hold on;
        plot(peakLocation(i)./FS,-projectedData(peakLocation(i)),'or','markerSize',10);
        hold off
        %axes(ax(3)); hold on; line([peakLocation(1)./500;peakLocation(1)./500],0.05,'linewidth',1.5, 'color','r');
        %hold off;
        %pause(0.2)
        ccc=input('Press enter: ','s');
        if ccc=='c'
            GRAPH=0;
        end
    end
    figure(fig);
    drawnow
end

%-------------------------------------------------------------------------%
%--------------------------------Update video-----------------------------%
%-------------------------------------------------------------------------%

if VIDEO
    %Coverage
    set(el(1:8),'FaceColor',[0 0 0]);
    set(el(sensor_coverage==1),'FaceColor',[0 0.7 0.1]);
    set(el(sensor_coverage==2),'FaceColor',[0.5 0.7 0]);
    set(el(sensor_coverage==3),'FaceColor',[1 0 0]);
    
    colormap([gray(64);jet(64)])
    
    %Camera
    imCAM = mcam.data(firstFrameCAM).bdata;
    imCAM=imCAM';
    imCAM=imCAM(end:-1:1,:); %for trail 2 only
    imCAM=imCAM(:,end:-1:1);
    imCAM=double(imCAM);
    imCAM=imCAM/max(imCAM(:))*63+1;
    image(imCAM,'Parent',ha(1),'CDataMapping','scaled'); %for trial 1, no rotation
    set(ha(1),'Xtick',[],'Ytick',[],'DataAspectRatio',[1 1 1],'CLim',[1 128]);
    
    %Pressure
    imPRESS = mpress.data(firstFramePRESS).bdata; %should be implemented in segments!
    imPRESS=imPRESS';
    imPRESS=imPRESS(end:-1:1,:); %for trial 1 and probably 2
    imPRESS=double(imPRESS);
    imPRESS=(imPRESS/max(imPRESS(:))*63)+65;
    a2=image(imPRESS,'Parent',ha(2),'CDataMapping','scaled');
    %surf(ha(2),imPRESS);
    set(ha(2),'Xtick',[],'Ytick',[],'DataAspectRatio',[1 1 1],'CLim',[1 128]);
    
    % Display ECG
    plot(ha(3),time,filteredDataREF, 'r');
    ylim(ha(3),[-400000 600000]); %-4 6, -5 8
    if RESP
        plot(ha(7),time,dataResp,'lineWidth',1.7,'color',[0 0.75 0.75]);
        %ylim(ha(7),[42100 43300]);
    end
    if KALMAN
        plot(ha(4),time,-cleanProjData,'b');
    else
        plot(ha(4),time,-projectedData),'b';
    end
    ylim(ha(4),[-10000 10000]); %-1 1 , -4 4
    
    % Display heart rate
    if numel(RRpositions)>=2
        irr=diff(RRpositions);
        startpos=max(2,find(RRpositions>=RRpositions(end-1)-15000,1));
        xlim('manual')
        plot(ha(6), RRpositions(startpos:end)./FS,60*FS./irr(startpos-1:end),'.b');
        if RRpositions(startpos)~=RRpositions(end)
            xlim(ha(6),[RRpositions(startpos)./FS, RRpositions(end)./FS]);
        end
    end
    %         if numel(RRpositionsALT)>=2  %AS
    %             axes(ha(6));
    %             hold on
    %             irr=diff(RRpositionsALT);
    %             plot(ha(6), RRpositionsALT(max(2,find(RRpositionsALT>=RRpositionsALT(end-1)-15000)):end)./FS,60*FS./irr(max(2,find(RRpositionsALT>=RRpositionsALT(end-1)-15000))-1:end),'.g');
    %         end
    if numel(RRpositionsREF)>=2
        %axes(ha(5));
        %hold on
        irr=diff(RRpositionsREF);
        startpos=max(2,find(RRpositionsREF>=RRpositionsREF(end-1)-15000,1));
        plot(ha(5), RRpositionsREF(startpos:end)./FS,60*FS./irr(startpos-1:end),'.r');
        if RRpositionsREF(startpos)~=RRpositionsREF(end)
            xlim(ha(5),[RRpositionsREF(startpos)./FS, RRpositionsREF(end)./FS]);
        end
        %hold off;
    end
    ylim(ha(5),[115 185]); %video martijn 100 150, 130 180, 110 160
    ylim(ha(6),[115 185]);
    set(ha(3),'Xtick',[],'Ytick',[]);
    set(ha(4),'Xtick',[],'Ytick',[]);
    set(ha(5),'Xtick',[]);
    set(ha(6),'Xtick',[]);
    title(ha(2),'Pressure distribution',...
        'FontName','Tw Cen MT',...
        'FontSize',15);
    title(ha(3),'Reference ECG',...
        'FontName','Tw Cen MT',...
        'FontSize',15);
    if RESP
        set(ha(7),'Xtick',[],'Ytick',[]);
        title(ha(7),'Contactless Respiration',...
            'FontName','Tw Cen MT',...
            'FontSize',15);
    end
    title(ha(4),'Contactless ECG',...
        'FontName','Tw Cen MT',...
        'FontSize',15);
    title(ha(5),'Heart Rate',...
        'FontName','Tw Cen MT',...
        'FontSize',15);
    % if sampleCounter >= 65*8000 && sampleCounter <= 95*8000 %video martijn chest
    % if sampleCounter >=   2064000 && sampleCounter <= 2400000 %back no data
    % if sampleCounter >= 3664000  && sampleCounter <=  3936000 %back
    writeVideo(vidObj,getframe(figVid));
    % end
end

%-------------------------------------------------------------------------%
%------------------------Get motion level based on video------------------%
%-------------------------------------------------------------------------%

if CAMERA
    imCAM = mcam.data(firstFrameCAM).bdata;
    motionLevel=MotionEstimator(imCAM);
    if SAVE
        safeMotionLevel=[safeMotionLevel, motionLevel];
    end
end

%-------------------------------------------------------------------------%
%------------------------Get motion level based on pressure mat-----------%
%-------------------------------------------------------------------------%

if PRESSMAT
    %imPRESS=cell(1,40);
    motionLevel=zeros(1,(40/8000*WINDOW_SIZE));
    % Get 40 frames per second (20 per half second) and take the squared diff
    % between 2 frames
    for i=1:(40/8000*WINDOW_SIZE)
        diff_frame = (mpress.data(min(firstFramePRESS+i,sizeFilePRESS)).bdata-mpress.data(firstFramePRESS+i-1).bdata).^2;
        motionLevel(i) = sum(sum(diff_frame));
    end
    if SAVE
        safeMotionLevel=[safeMotionLevel, var(motionLevel)];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%             Remarks about the data processing                %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Search for all the "AS" to find related comments...
