% LAUNCH_PROCESSING.M
%
% This file:
% - defines the conditions in which the main code should run for offline data analysis;
% - defines the dataset to be analyzed;
% - is not needed for the online demo on LabView;
%
% This file:
% - is called directly from the MATLAB workspace;
% - calls the file read_data.m
%
% aline.serteyn@gmail.com
% August 2012
%
% -------------------------------------------
% For online processing of online data (DEMO):
% -------------------------------------------
%
% - start LabView and the E-NEMO project
%
% - make the LabView/Matlab interface call main_algorithm when processing is needed
%
% - run launch_online_demo in the Matlab command window
%
% - start the demo
%
%
% -------------------------------------
% For online processing of offline data:
% -------------------------------------
%
% - start Matlab
%
% - make sure the current directory is the cECG toolbox
%
% (- Optional: define the processing mode in launch_processing.m (video yes or not, dynamic graph yes or not,specific dataset...))
%
% (- Optional: edit the file launch_processing.m to further process the data or display specific graphs)
%
% - run launch_processing

% clear;
tic


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                Add the required directories                  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add processing toolboxes
% ------------------------
CurrentDirectory = cd;
addpath(cd)
addpath([cd '\cECG_processing_functions'])
addpath([cd '\cECG_display_functions'])


% Localize and group datasets
% ---------------------------

% >>>>>>>>>> Baby 4 on his chest <<<<<<<<<<<<<<<<
% ADULT=0;
% Body=1;
% dataSet{1} = 'D:\Philips DATA\20120704_test4\20120704_143614'; numFiles_dataSet(1)=2; STARTtime_dataSet(1)=45; ENDtime_dataSet(1)=49;
% dataSet{2} = 'D:\Philips DATA\20120704_test4\20120704_143614'; numFiles_dataSet(2)=3; STARTtime_dataSet(2)=10; ENDtime_dataSet(2)=26;

% >>>>>>>>>> Critical periods <<<<<<<<<<<<<<<<
% on chest:
ADULT=0;
Body=1;
%The follwoing directory is where you're saving processed files; This is
%determined in C_create_cECG_mov located in ... Matlab\cECG Data specific
savingdir=dataSetpath{3};
dataSet{1} = dataPath; numFiles_dataSet(1)= numFiles; STARTtime_dataSet(1)=0; ENDtime_dataSet(1)=Inf; savedir{1}= outputPath;

totalTime_iHR_WD=0; %in minutes
goodCoverTime_iHR_WD=0;
% totalTime_iHR_WD2=0;
% goodCoverTime_iHR_WD=0;
totalTime_iHR_NI=0;
totalTime_iHR_NM=0;
totalTime_iHR_LM=0;
totalTime_iHR_HM=0;
totalTime_iHR_IN=0;
goodCoverTime_iHR_NI=0;
numSeg_iHR_NI=0;
maxSeg_iHR_NI=0;
ppv_iHR_WD=0;
sen_iHR_WD=0;
err_iHR_WD=0;
ppv_iHR_NI=0;
sen_iHR_NI=0;
err_iHR_NI=0;
ppv_iHR_NM=0;
sen_iHR_NM=0;
err_iHR_NM=0;
ppv_iHR_LM=0;
sen_iHR_LM=0;
err_iHR_LM=0;
ppv_iHR_HM=0;
sen_iHR_HM=0;
err_iHR_HM=0;
ppv_iHR_IN=0;
sen_iHR_IN=0;
err_iHR_IN=0;
%totalTime_iHR_NI2=0;
goodCoverTime_iHR_NI2=0;
totalTime_aHR_NI0=0;
goodTime_aHR_NI=0;
badTime_aHR_NI=0;
goodTime_aHR_NI0=0;
totalTime_aHR_WD0=0;
goodTime_aHR_WD=0;
badTime_aHR_WD=0;
goodTime_aHR_WD0=0;
totalTime_aHR_NM=0;
goodTime_aHR_NM=0;
badTime_aHR_NM=0;
totalTime_aHR_LM=0;
goodTime_aHR_LM=0;
badTime_aHR_LM=0;
totalTime_aHR_HM=0;
goodTime_aHR_HM=0;
badTime_aHR_HM=0;
totalTime_aHR_IN=0;
goodTime_aHR_IN=0;
badTime_aHR_IN=0;
numSeg_aHR_NI=0;
maxSeg_aHR_NI=0;
numSeg_aHR_NI0=0;
maxSeg_aHR_NI0=0;
non_matching_segments=[];

for index_dataset=1:numel(dataSet)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%                Define the processing mode (booleans)         %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Fixed values:
    % -------------
    ONLINE_DEMO=0; %=1 for the LabView online demo --> ALWAYS 1 here
    DO_INIT=1; %ALWAYS 1 here: initialize other variables
    KALMAN_RESET=1; %ALWAYS 1 here:initialize kalman variables
    PEAK_RESET=1; %ALWAYS 1 here:initialize peak detection variables
    Manual=1; %body position determination (always 1 until the automated one is implemented)
    
    % May be changed:
    % --------------
    % default: 4000, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1
    WINDOW_SIZE=8000;%8000; % in samples; min 352 (delay of pre-resampling filter) and max 16000 (buffer size) for 8000Hz, (Note: RRmin = 250bpm = 1900 samples @ 8000Hz)
    VIDEO=0; %creates a video
    RESP=0; %compute the respiration
    KALMAN=0; %aply Kalman smoother to the data
    GRAPH=0; %display online graphs => should be 0 for online demo
    CAMERA=0; %read the video files for motion level estimation
    if CAMERA; addpath('C:\Users\Aline\Documents\Research\Philips cECG\MATLAB Philips\Marek\motion_detection');end
    PRESSMAT=0; %read the pressure mat files for motion level estimation
    ANNOTATIONS=0; %read the annotations (.xlsx files) for motion level estimation
    VERBOSE=1; %display comment and info about the progress of the algorithm
    SAVE=1; %concatenates processed windows (=get full processed signals) for later use => should be 0 for online demo to prevent memory overflow
    FIXED_NUM_RR=0; %choice of averaging strategy for the average heart rate (meanRR): 1=fixed number of RRintervals used for the computation of the meanRR, 0=last 10seconds (at least) of data used for the computation of the meanRR
    ALTERNATIVE_PEAK_DETECT=1; %if ALTERNATIVE_PEAK_DETECT=1 --> detect peaks on an alternative signal (e.g. diff of 2 fixed sensors) in parallel with the normal processing (used for offline analysis)
    CHAN_SELECT=3; %CHAN_SELECT=0=old Aline selection strategy, CHAN_SELECT=1=new Louis selection strategy, CHAN_SELECT=[]=all 8 channels selected, CHAN_SELECT=3 = 2channels selected based on corr with REF
    CHAN_2SELECT=1; %WHEN ALTERNATIVE_PEAK_DETECT=1 --> CHAN_2SELECT=0=always #2-#4; CHAN_2SELECT=1=always 2 best coupled channels (aline''s or Louis); CHAN_2SELECT=2=always 2 best non-adjacents
    MIN3CHANNELS=0; %force the selection of minimum 3 covered (not necessarly nicely coupled) channels in order to compute the 3 Einthoven leads
    NOTCH24=1; %if =1: applies a 24Hz notch filter to the data
    SAVE_TO_DIR = 1; %save the cap ecg and HR to a given directory. This is done per dataset, hence per recording session (per data folder).
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%                Launch processing                             %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('-------------------------------------------------------------------')
    disp(['Processing data set:   ' dataSet{index_dataset}])
    disp('-------------------------------------------------------------------')
    DataDirectory=dataSet{index_dataset}; numFiles=numFiles_dataSet(index_dataset); STARTtime=STARTtime_dataSet(index_dataset); ENDtime=ENDtime_dataSet(index_dataset);
    % *** edit by jan 
    cd('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Create data\creating cECG and movement files\cECG_processing_functions')
    % *** edit by jan 
    
    
% *****************************************
% *****************************************
if jansway==0
    read_data % original
elseif jansway==1
    read_data_jansway % changed by Jan. Should fix the multiple file issue
end
% *****************************************
% *****************************************

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%                Update global results                         %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Note: check the structures results and resultsALT!
    % results is for the CHAN_SELECT strategy
    % resultsALT is for the CHAN_SELECT2 strategy
    if SAVE
        % coverage for the whole data with several channels selected
        totalTime_iHR_WD=totalTime_iHR_WD+results.wholeData(1); %in minutes
        goodCoverTime_iHR_WD=goodCoverTime_iHR_WD+results.wholeData(2); %in minutes
        
        % total time of the whole data (intervention and ignored data excluded) with several channels selected
        totalTime_iHR_NI=totalTime_iHR_NI+results.wholeData(1)-results.intervention(1)-results.ignored(1); %in minutes
        % coverage for the whole data (intervention and ignored data excluded) with several channels selected
        goodCoverTime_iHR_NI=goodCoverTime_iHR_NI+results.wholeData(2)-results.intervention(2)-results.ignored(2); %in minutes
        % number of segments > 10s for the whole data (intervention and ignored data excluded) with several channels selected
        numSeg_iHR_NI=numSeg_iHR_NI+results.noMotion(3)+results.lowMotion(3)+results.highMotion(3)+results.ghosts(3);
        maxSeg_iHR_NI=max([maxSeg_iHR_NI, results.noMotion(4), results.lowMotion(4), results.highMotion(4), results.ghosts(4)]);
        % ppv for the whole data (intervention and ignored data excluded) with several channels selected
        ppv_iHR_NI=ppv_iHR_NI+results.wholeData(1)*results.wholeData(6)-results.intervention(1)*results.intervention(6)-results.ignored(1)*results.ignored(6); %in minutes
        % sensitivity for the whole data (intervention and ignored data excluded) with several channels selected
        sen_iHR_NI=sen_iHR_NI+results.wholeData(1)*results.wholeData(7)-results.intervention(1)*results.intervention(7)-results.ignored(1)*results.ignored(7); %in minutes
        % error for the whole data (intervention and ignored data excluded) with several channels selected
        err_iHR_NI=err_iHR_NI+results.wholeData(1)*results.wholeData(8)-results.intervention(1)*results.intervention(8)-results.ignored(1)*results.ignored(8); %in minutes
        % Whole data
        ppv_iHR_WD=ppv_iHR_WD+results.wholeData(1)*results.wholeData(6); %in minutes
        sen_iHR_WD=sen_iHR_WD+results.wholeData(1)*results.wholeData(7); %in minutes
        err_iHR_WD=err_iHR_WD+results.wholeData(1)*results.wholeData(8); %in minutes
        % No motion
        totalTime_iHR_NM=totalTime_iHR_NM+results.noMotion(1);
        ppv_iHR_NM=ppv_iHR_NM+results.noMotion(1)*results.noMotion(6); %in minutes
        sen_iHR_NM=sen_iHR_NM+results.noMotion(1)*results.noMotion(7); %in minutes
        err_iHR_NM=err_iHR_NM+results.noMotion(1)*results.noMotion(8); %in minutes
        % Low motion
        totalTime_iHR_LM=totalTime_iHR_LM+results.lowMotion(1);
        ppv_iHR_LM=ppv_iHR_LM+results.lowMotion(1)*results.lowMotion(6); %in minutes
        sen_iHR_LM=sen_iHR_LM+results.lowMotion(1)*results.lowMotion(7); %in minutes
        err_iHR_LM=err_iHR_LM+results.lowMotion(1)*results.lowMotion(8); %in minutes
        % High motion
        totalTime_iHR_HM=totalTime_iHR_HM+results.highMotion(1);
        ppv_iHR_HM=ppv_iHR_HM+results.highMotion(1)*results.highMotion(6); %in minutes
        sen_iHR_HM=sen_iHR_HM+results.highMotion(1)*results.highMotion(7); %in minutes
        err_iHR_HM=err_iHR_HM+results.highMotion(1)*results.highMotion(8); %in minutes
        % Intervention
        totalTime_iHR_IN=totalTime_iHR_IN+results.intervention(1);
        ppv_iHR_IN=ppv_iHR_IN+results.intervention(1)*results.intervention(6); %in minutes
        sen_iHR_IN=sen_iHR_IN+results.intervention(1)*results.intervention(7); %in minutes
        err_iHR_IN=err_iHR_IN+results.intervention(1)*results.intervention(8); %in minutes
        
        
        % HISTOGRAM
        non_matching_segments=[non_matching_segments, results.wholeData(9:end)];
        
        
        if ALTERNATIVE_PEAK_DETECT
            %     % coverage for the whole data with only 2 channels selected
            %     totalTime_iHR_WD2=totalTime_iHR_WD2+resultsALT.wholeData(1);
            %     goodCoverTime_iHR_WD2=goodCoverTime_iHR_WD2+resultsALT.wholeData(2);
            % coverage for the whole data (intervention and ignored data excluded) with only 2 channels selected
            %totalTime_iHR_NI2=totalTime_iHR_NI2+resultsALT.wholeData(1)-resultsALT.intervention(1)-resultsALT.ignored(1); %in minutes
            goodCoverTime_iHR_NI2=goodCoverTime_iHR_NI2+resultsALT.wholeData(2)-resultsALT.intervention(2)-resultsALT.ignored(2); %in minutes
        end
        % Max segment/ mean segment/ num of segments (3bpm)
        
        %>> average HR coverage (3bpm)
        %%% without rel indicator:
        % time coverage:
        % ...whole data...
        error_temp=error0;
        totalTime_aHR_WD0=totalTime_aHR_WD0+length(error_temp)/8000*WINDOW_SIZE/60; %in minutes
        goodTime_aHR_WD0=goodTime_aHR_WD0+sum(error_temp<=3)/8000*WINDOW_SIZE/60; %in minutes
        % ...no intervention...
        error_temp=error0(motionIn==1|motionIn==2|motionIn==3|motionIn==6);
        totalTime_aHR_NI0=totalTime_aHR_NI0+length(error_temp)/8000*WINDOW_SIZE/60; %in minutes
        goodTime_aHR_NI0=goodTime_aHR_NI0+sum(error_temp<=3)/8000*WINDOW_SIZE/60; %in minutes
        % segment length when no intervention:
        segLength0=evaluate_matching([], error_temp, 3, WINDOW_SIZE); %note that the error signal is cut so periods may appear shorter that they are in reality
        numSeg_aHR_NI0=numSeg_aHR_NI0+numel(segLength0);
        maxSeg_aHR_NI0=max([maxSeg_aHR_NI0,segLength0]);
        
        %%% with reliability indicator
        % time coverage
        % ...whole data...
        error_temp=errorRelInd;
        goodTime_aHR_WD=goodTime_aHR_WD+sum(error_temp<=3)/8000*WINDOW_SIZE/60;
        badTime_aHR_WD=badTime_aHR_WD+sum(error_temp<=Inf)/8000*WINDOW_SIZE/60;
        % ...no motion...
        error_temp=errorRelInd(motionIn==1);
        goodTime_aHR_NM=goodTime_aHR_NM+sum(error_temp<=3)/8000*WINDOW_SIZE/60;%in minutes
        badTime_aHR_NM=badTime_aHR_NM+sum(error_temp<=Inf)/8000*WINDOW_SIZE/60;%in minutes
        totalTime_aHR_NM=totalTime_aHR_NM+length(error_temp)/8000*WINDOW_SIZE/60; %in minutes
        % ...low motion...
        error_temp=errorRelInd(motionIn==2);
        goodTime_aHR_LM=goodTime_aHR_LM+sum(error_temp<=3)/8000*WINDOW_SIZE/60;%in minutes
        badTime_aHR_LM=badTime_aHR_LM+sum(error_temp<=Inf)/8000*WINDOW_SIZE/60;%in minutes
        totalTime_aHR_LM=totalTime_aHR_LM+length(error_temp)/8000*WINDOW_SIZE/60; %in minutes
        % ...high motion...
        error_temp=errorRelInd(motionIn==3);
        goodTime_aHR_HM=goodTime_aHR_HM+sum(error_temp<=3)/8000*WINDOW_SIZE/60;%in minutes
        badTime_aHR_HM=badTime_aHR_HM+sum(error_temp<=Inf)/8000*WINDOW_SIZE/60;%in minutes
        totalTime_aHR_HM=totalTime_aHR_HM+length(error_temp)/8000*WINDOW_SIZE/60; %in minutes
        % ...intervention...
        error_temp=errorRelInd(motionIn==4|motionIn==5);
        goodTime_aHR_IN=goodTime_aHR_IN+sum(error_temp<=3)/8000*WINDOW_SIZE/60;%in minutes
        badTime_aHR_IN=badTime_aHR_IN+sum(error_temp<=Inf)/8000*WINDOW_SIZE/60;%in minutes
        totalTime_aHR_IN=totalTime_aHR_IN+length(error_temp)/8000*WINDOW_SIZE/60; %in minutes
        % ...no intervention...
        error_temp=errorRelInd(motionIn==1|motionIn==2|motionIn==3|motionIn==6);
        goodTime_aHR_NI=goodTime_aHR_NI+sum(error_temp<=3)/8000*WINDOW_SIZE/60;%in minutes
        badTime_aHR_NI=badTime_aHR_NI+sum(error_temp<=Inf)/8000*WINDOW_SIZE/60;%in minutes
        % segment length when no intervention
        segLength=evaluate_matching([], error_temp, 3, WINDOW_SIZE); %note that the error signal is cut so periods may appear shorter that they are in reality
        numSeg_aHR_NI=numSeg_aHR_NI+numel(segLength);
        maxSeg_aHR_NI=max([maxSeg_aHR_NI,segLength]);
        
        
        %Blant-Altman plot for this data directory
        %aCTIVATED BY jAN 
        %     %Saving matrices
            matrixA{index_dataset}=[timeWin',position',motionIn',(safeMeanRR')./FS,(safeMeanRR2')./FS,(safeMeanRRREF')./FS]; %var_ampli2 %RR intervals in seconds
            all_cleanProjECG{index_dataset}=cleanProjectedECG;
            all_projECG{index_dataset}=projectedECG;
            all_RRpositions{index_dataset}=RRpositions;
            all_RRpositionsREF{index_dataset}=RRpositionsREF;
        
        if SAVE_TO_DIR
            motion_level = var_ampli2; %FS Hz
            save([savedir{index_dataset} '\motion_level.mat'], 'motion_level')
            cap_ECG = projectedECG; % FS Hz
            save([savedir{index_dataset} '\cap_ECG.mat'], 'cap_ECG')
            cap_Rpeakposition = RRpositions; %in samples
            save([savedir{index_dataset} '\cap_Rpeakposition.mat'], 'cap_Rpeakposition')
            ref_Rpeakposition = RRpositionsREF; %in samples
            save([savedir{index_dataset} '\ref_Rpeakposition.mat'], 'ref_Rpeakposition')
            save([savedir{index_dataset} '\referenceECG.mat'], 'referenceECG')
        end
        
        %save('matrixA.m','matrixA')
        %figure; plot(abs(safeMeanRRREF./FS-safeMeanRR./FS),mean([safeMeanRRREF./FS;safeMeanRR./FS],1),'.'); title('Bland-Altman plot')
        %NOTE: BA plot should be from RR tacogram not the averaged RR intervals!!! --> need to evenly sample the RR tacogram (or iHR signal)+ discretise with respect to motion (color code)!!
    end
end

cd(CurrentDirectory)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                Initialize text file                          %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SAVE && ~SAVE_TO_DIR
    filename = sprintf('Results_%s.txt', 'several_babies'); % number serves only as an example
    fid2 = fopen(filename,'a+t'); %w+t, open or create file (in text mo) for reading and writing; discard existing contents!
    fprintf(fid2,'%s\n','---------------------------------------------------------------------------------------------------------------------------');
    fprintf(fid2,'%s\t\t\t%s\n','Analysis date:', datestr(now,21)); % \t means "tab space"
    fprintf(fid2,'%s\t\t\t', 'Data set: ');
    for i=1:numel(dataSet)
        fprintf(fid2,'%s; ', [dataSet{i}(end-14:end) ' (' num2str(STARTtime_dataSet(i)) '-' num2str(ENDtime_dataSet(i)) ')']);
    end
    fprintf(fid2,'\n');
    fprintf(fid2,'%s\t\t%s\n','Sampling frequency:', num2str(FS));
    %fprintf(fid2,'%s\t%s%s\n','Length processed signal:', num2str(sampleCounter/8000/60-initOffset/60), ' minutes');
    %fprintf(fid2,'%s\t\t\t%s%s\n','Starting time:', num2str(initOffset/60), ' minutes after the beginning of the recording.');
    if ALTERNATIVE_PEAK_DETECT
        fprintf(fid2,'%s\t\t\n',['Configuration: WINDOW_SIZE=' num2str(WINDOW_SIZE)  ',KALMAN=' num2str(KALMAN) ',MIN3CHANNELS=' num2str(MIN3CHANNELS) ',CHAN_SELECT=' num2str(CHAN_SELECT) ',CHAN_2SELECT=' num2str(CHAN_2SELECT) ',NOTCH24=' num2str(NOTCH24) ',FIXED_NUM_RR=' num2str(FIXED_NUM_RR) ',ANNOTATIONS=' num2str(ANNOTATIONS) ',RESP=' num2str(RESP) ',PRESSMAT=' num2str(PRESSMAT)]);
    else
        fprintf(fid2,'%s\t\t\n',['Configuration: WINDOW_SIZE=' num2str(WINDOW_SIZE)  ',KALMAN=' num2str(KALMAN) ',MIN3CHANNELS=' num2str(MIN3CHANNELS)  ',CHAN_SELECT=' num2str(CHAN_SELECT) ',NOTCH24=' num2str(NOTCH24)  ',FIXED_NUM_RR=' num2str(FIXED_NUM_RR) ',ANNOTATIONS=' num2str(ANNOTATIONS) ',RESP=' num2str(RESP) ',PRESSMAT=' num2str(PRESSMAT)]);
    end
    if ANNOTATIONS
        fprintf(fid2,'%s\n','Motion segmentation based on: annotations');
    elseif RESP
        fprintf(fid2,'%s\n','Motion segmentation based on: 1 kHz (Louis)');
    elseif PRESSMAT
        fprintf(fid2,'%s\n','Motion segmentation based on: pressure-mat');
    else
        fprintf(fid2,'%s\n','Motion segmentation based on: 1 kHz (Aline)');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%                Save/display global results                   %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf(fid2,'%s\n','---------------------------------------------------------------------------------------------------------------------------');
    fprintf(fid2,'%s\n',['Results for the INSTANTANEOUS heart rate (peak to peak basis, with channel selection strategy #' num2str(CHAN_SELECT) ').']);
    fprintf(fid2,'%s\n','---------------------------------------------------------------------------------------------------------------------------');
    fprintf(fid2, '\n');
    % Time coverage
    fprintf(fid2, '%s\n','Time coverage:');
    fprintf(fid2, '\t%s\n',['total = ' num2str(totalTime_iHR_NI) ' min; good% = ' num2str(round(goodCoverTime_iHR_NI/totalTime_iHR_NI*100)) ' % (' num2str(round(goodCoverTime_iHR_NI2/totalTime_iHR_NI*100)) ' %)']);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s\n','Interpretation:');
    fprintf(fid2, '%s\n',['For ' num2str(totalTime_iHR_WD) ' minutes of data, the RR-intervals where perfectly detected during ' num2str(goodCoverTime_iHR_WD) ' minutes (i.e. ' num2str(round(goodCoverTime_iHR_WD/totalTime_iHR_WD*100)) '% of the time).']);
    fprintf(fid2, '%s\n',['If we exclude the interventions and ignored periods, we are left with ' num2str(totalTime_iHR_NI) ' minutes of data from which ' num2str(goodCoverTime_iHR_NI) ' minutes lead to a perfect detection of RR-intervals.']);
    fprintf(fid2, '%s\n',['In other words, when the baby doesn''t receive special care, the instantaneous heart rate is perfectly reliable during ' num2str(round(goodCoverTime_iHR_NI/totalTime_iHR_NI*100)) '% of the time.']);
    if ALTERNATIVE_PEAK_DETECT
        fprintf(fid2, '%s\n',['(Note that his percentage of time is changed to ' num2str(round(goodCoverTime_iHR_NI2/totalTime_iHR_NI*100)) '% if we use only 2 sensors from the array of 8 sensors.)']);
    end
    % Segments > 10s
    fprintf(fid2, '\n');
    fprintf(fid2, '%s\n','Lost segments:');
    fprintf(fid2, '\t%s\n',['# = ' num2str(numSeg_iHR_NI) ' ; max = ' num2str(maxSeg_iHR_NI) ' s']);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s\n','Interpretation:');
    fprintf(fid2, '%s\n', ['If we exclude the interventions and ignored periods, there were ' num2str(numSeg_iHR_NI) ' segments of data LONGER than 10 seconds which didn''t allowed an accurate instantaneous heart rate computation.']);
    fprintf(fid2, '%s\n', ['The longest period of unreliable instantaneous heart rate for this dataset (intervention and ignored periods excluded) was of ' num2str(maxSeg_iHR_NI) ' seconds.']);
    fprintf(fid2, '\n');
    % Details per motion
    fprintf(fid2, '%s\n', 'Details of the peak detection for different motion levels:');
    fprintf(fid2,'\t%s\t%s\t%s\t%s\t%s\n', 'Activity', 'Minutes', 'PPV', 'SEN', 'ERROR');
    fprintf(fid2,'\t%s\n', '----------------------------------');
    fprintf(fid2,'\t%s\t%2.1f\t%2.1f\t%2.1f\t%1.2f\n', 'No motion', totalTime_iHR_NM, ppv_iHR_NM/totalTime_iHR_NM, sen_iHR_NM/totalTime_iHR_NM, err_iHR_NM/totalTime_iHR_NM);
    fprintf(fid2,'\t%s\t%2.1f\t%2.1f\t%2.1f\t%1.2f\n', 'Low motion', totalTime_iHR_LM, ppv_iHR_LM/totalTime_iHR_LM, sen_iHR_LM/totalTime_iHR_LM, err_iHR_LM/totalTime_iHR_LM);
    fprintf(fid2,'\t%s\t%2.1f\t%2.1f\t%2.1f\t%1.2f\n', 'High motion', totalTime_iHR_HM, ppv_iHR_HM/totalTime_iHR_HM, sen_iHR_HM/totalTime_iHR_HM, err_iHR_HM/totalTime_iHR_HM);
    fprintf(fid2,'\t%s\t%2.1f\t%2.1f\t%2.1f\t%1.2f\n', 'Intervent.', totalTime_iHR_IN, ppv_iHR_IN/totalTime_iHR_IN, sen_iHR_IN/totalTime_iHR_IN, err_iHR_IN/totalTime_iHR_IN);
    fprintf(fid2,'\t%s\n', '----------------------------------');
    fprintf(fid2,'\t%s\t%2.1f\t%2.1f\t%2.1f\t%1.2f\n', 'Whole data', totalTime_iHR_WD, ppv_iHR_WD/totalTime_iHR_WD, sen_iHR_WD/totalTime_iHR_WD, err_iHR_WD/totalTime_iHR_WD);
    fprintf(fid2,'\t%s\n', '----------------------------------');
    fprintf(fid2,'\t%s\t%2.1f\t%2.1f\t%2.1f\t%1.2f\n', 'No intervent.', totalTime_iHR_NI, ppv_iHR_NI/totalTime_iHR_NI, sen_iHR_NI/totalTime_iHR_NI, err_iHR_NI/totalTime_iHR_NI);
    fprintf(fid2, '\n');
    
    
    if ~FIXED_NUM_RR
        fprintf(fid2,'%s\n','---------------------------------------------------------------------------------------------------------------------------');
        fprintf(fid2,'%s\n',['Results for the AVERAGED heart rate (mean iHR value over the last ' num2str(NUM_RR_AVERAGED) ' seconds (at least)).']);
        fprintf(fid2,'%s\n','---------------------------------------------------------------------------------------------------------------------------');
    else
        fprintf(fid2,'%s\n','---------------------------------------------------------------------------------------------------------------------------');
        fprintf(fid2,'%s\n',['Results for the AVERAGED heart rate (mean iHR value over the last ' num2str(NUM_RR_AVERAGED) ' RR intervals).']);
        fprintf(fid2,'%s\n','---------------------------------------------------------------------------------------------------------------------------');
    end
    % Time coverage
    fprintf(fid2, '\n');
    fprintf(fid2, '%s\n','Time coverage:');
    fprintf(fid2, '\t%s\n',['total = ' num2str(totalTime_aHR_NI0) ' min; good% = ' num2str(round(goodTime_aHR_NI0/totalTime_aHR_NI0*100)) ' %, bad% = ' num2str(100-round(goodTime_aHR_NI0/totalTime_aHR_NI0*100))]);
    fprintf(fid2, '\t%s\n',['Rel. Indic. ----> good% = ' num2str(round(goodTime_aHR_NI/totalTime_aHR_NI0*100)) ' %, bad% = ' num2str(round((badTime_aHR_NI-goodTime_aHR_NI)/totalTime_aHR_NI0*100))]);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s\n','Interpretation:');
    fprintf(fid2, '%s\n',['By comparing the averaged HR signal with the reference averaged HR signal, an error may occur.']);
    fprintf(fid2, '%s\n', ['If we assume that an error of 3 bpm is acceptable, then we have an reliable heart rate on the whole dataset for ' num2str(round(goodTime_aHR_WD0/totalTime_aHR_WD0*100)) '% of the ' num2str(totalTime_aHR_WD0) ' minutes considered.']);
    fprintf(fid2, '%s\n',['(Note that if a reliability indicator is used during the computation of the mean, then we get a reliable heart rate on the whole dataset for ' num2str(round(goodTime_aHR_WD/totalTime_aHR_WD0*100)) '% of the ' num2str(totalTime_aHR_WD0) ' minutes considered.)']);
    fprintf(fid2, '%s\n', ['If we assume that an error of 3 bpm is acceptable and if we exclude the intervention periods, then we have an reliable heart rate for ' num2str(round(goodTime_aHR_NI0/totalTime_aHR_NI0*100)) '% of the ' num2str(totalTime_aHR_NI0) ' minutes considered.']);
    fprintf(fid2, '%s\n',['(Note that if a reliability indicator is used during the computation of the mean, then we get a reliable heart rate (intervention excluded) for ' num2str(round(goodTime_aHR_NI/totalTime_aHR_NI0*100)) '% of the ' num2str(totalTime_aHR_NI0) ' minutes considered.)']);
    % Segments > 10s
    fprintf(fid2, '\n');
    fprintf(fid2, '%s\n','Lost segments:');
    fprintf(fid2, '\t%s\n',['# = ' num2str(numSeg_aHR_NI0) '; max = ' num2str(maxSeg_aHR_NI0) ' s']);
    fprintf(fid2, '\t%s\n',['Rel. Indic. --> # = ' num2str(numSeg_aHR_NI) '; max = ' num2str(maxSeg_aHR_NI) ' s']);
    fprintf(fid2, '\n');
    fprintf(fid2, '%s\n','Interpretation:');
    fprintf(fid2, '%s\n', ['If we exclude the interventions and ignored periods, there were ' num2str(numSeg_aHR_NI0) ' segments of data LONGER than 10 seconds during which the averaged heart rate was not matching the reference averaged heart rate.']);
    fprintf(fid2, '%s\n', ['The longest period of unreliable averaged heart rate for this dataset (intervention and ignored periods excluded) was of ' num2str(maxSeg_aHR_NI0) ' seconds.']);
    fprintf(fid2, '%s\n', ['(Note that if a reliability indicator is used during the computation of the mean, we got ' num2str(numSeg_aHR_NI) ' segments of data LONGER than 10 seconds whose max is of ' num2str(maxSeg_aHR_NI) ' seconds.)']);
    fprintf(fid2, '\n');
    % Per motion
    fprintf(fid2, '%s\n', 'Details of the averaged heart rate for different motion levels (with Rel. Ind. and 3bpm difference allowed):');
    fprintf(fid2,'\t%s\t\t%s\t\t%s\t\t%s\n', 'Activity', 'Minutes', 'good%', 'bad%');
    fprintf(fid2,'\t%s\n', '-------------------------------------------------------');
    fprintf(fid2,'\t%s\t\t%3.2f\t\t%2i\t\t%2i\n', 'No motion   ', totalTime_aHR_NM, round(goodTime_aHR_NM/totalTime_aHR_NM*100), round((badTime_aHR_NM-goodTime_aHR_NM)/totalTime_aHR_NM*100));
    fprintf(fid2,'\t%s\t\t%3.2f\t\t%2i\t\t%2i\n', 'Low motion  ', totalTime_aHR_LM, round(goodTime_aHR_LM/totalTime_aHR_LM*100), round((badTime_aHR_LM-goodTime_aHR_LM)/totalTime_aHR_LM*100));
    fprintf(fid2,'\t%s\t\t%3.2f\t\t%2i\t\t%2i\n', 'High motion ', totalTime_aHR_HM, round(goodTime_aHR_HM/totalTime_aHR_HM*100), round((badTime_aHR_HM-goodTime_aHR_HM)/totalTime_aHR_HM*100));
    fprintf(fid2,'\t%s\t\t%3.2f\t\t%2i\t\t%2i\n', 'Intervention', totalTime_aHR_IN, round(goodTime_aHR_IN/totalTime_aHR_IN*100), round((badTime_aHR_IN-goodTime_aHR_IN)/totalTime_aHR_IN*100));
    % error_temp=errorRelInd(motionIn==6);
    % fprintf(fid,'%s\t\t%3.2f\t\t\t%2i\t\t%2i\t\t%2i\n', 'Ghosts      ', sum(motionIn==6)/8000*WINDOW_SIZE/60, round(sum(error_temp<=1)/length(error_temp)*100), round(sum(error_temp<=3)/length(error_temp)*100), round(sum(error_temp<=5)/length(error_temp)*100));
    % error_temp=errorRelInd(motionIn==7);
    % fprintf(fid,'%s\t\t%3.2f\t\t\t%2i\t\t%2i\t\t%2i\n', 'Ignored Data', sum(motionIn==7)/8000*WINDOW_SIZE/60, round(sum(error_temp<=1)/length(error_temp)*100), round(sum(error_temp<=3)/length(error_temp)*100), round(sum(error_temp<=5)/length(error_temp)*100));
    fprintf(fid2,'\t%s\n', '-------------------------------------------------------');
    error_temp=errorRelInd;
    fprintf(fid2,'\t%s\t\t%3.2f\t\t%2i\t\t%2i\n', 'Whole data  ', totalTime_aHR_WD0, round(goodTime_aHR_WD/totalTime_aHR_WD0*100), round((badTime_aHR_WD-goodTime_aHR_WD)/totalTime_aHR_WD0*100));
    fprintf(fid2,'\t%s\n', '-------------------------------------------------------');
    fprintf(fid2,'\t%s\t\t%3.2f\t\t%2i\t\t%2i\n', 'No intervention', totalTime_aHR_NI0, round(goodTime_aHR_NI/totalTime_aHR_NI0*100), round((badTime_aHR_NI-goodTime_aHR_NI)/totalTime_aHR_NI0*100));
    fprintf(fid2, '\n');
    fprintf(fid2,'%s\n','---------------------------------------------------------------------------------------------------------------------------');
    
    fclose(fid2);
    
    
    % HISTOGRAM
    figure('Name', 'Histogram (iHR)');
    hist(non_matching_segments); title('Segment lengths distribution for iHR (whole data)');
    xlabel('Segment length (s)')
    ylabel('Number of segments')
    
end
cd(savingdir)
save matrixA matrixA;
save all_cleanProjECG all_cleanProjECG;
save all_projECG all_projECG;
save all_RRpositions all_RRpositions;
save all_RRpositionsREF all_RRpositionsREF;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%             Remarks about the data processing                %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% - Feel free to change the booleans VIDEO, RESP, GRAPH, CAMERA, VERBOSE and SAVE; as well as the processing window size WINDOW_SIZE

% - the minimal value for <numberNewSamples> is limited by min 50 (number coefficents of temporal filter) and min 80 (min number of samples provided by the recording system) and 352 (@8000Kz): delay of the preresampling filter and 160 for the search of the max in 1kHz

% - max 9 different binary files for ECG and RESP, 99 for pressure mat and 999 for video can be read synchronously using this interface

% - only one video frame is read per analysis window (WINDOW_SIZE) => if the VIDEO or CAMERA options are used, consider taking a small WINDOW_SIZE

% - the referenceECG of the new Philips monitors (from TRIAL2 onwards) is smoothed (it looks artificial) => it looks bad after bandpass filtering (ripples)

% - Impossible to create a video (read the camera and pressure mat data) for ECG data within the second ECG file! (synchronization problems)

% NOTE----->: the video created for Martijn's DOVO was from this dataset:
% DataDirectory = ' ...philips_share... \Trial 2\20120425_163034'; numFiles=1;ADULT=0; Body=1; STARTtime=1; ENDtime=2; VIDEO=1; RESP=1;
% <---------


% V28-08-2012:
% - Kalman resets only kalman variables!
% - MIN3CHANNELS: boolean
% - CHAN_2SELECT: for the peak detection strategy
% - FUSION renamed
% - text file improved with config

%V8-10-2012 (main, launch, determineFilter, notchFiltering, spotdiff,outpikeloc,flip_proximity):
% - Louis CHAN_SELECT strategy added (+CHAN_2SELECT)
% - Proximity_matrix init
% - notch filter init
% - debug CHAN_2SELECT==2: channelToUse instead of channelsToUseALT
% - notch filter (determineFilter + notchFiltering.m + filteredData1/filteredDataREF1)
% - move down safe block
% - vcg computations if num_channels>2
% - ATTENTION: rel_indicator: numel(channelsToUse)>=3 replaced by
% numel(channelsToUse)>=2 (line 840
% - histogram (offline proc, launch proc, determine_...)

%V6-04-2017
% added SAVE_TO_DIR to save the capacitive ECG for other analysis (e.g.
% slee stages with Jan Werth)
% NOTE: the delay between the reference ECG and the capacitive ECG coming
% out of the processing chain seems to be of 1.736 seconds (this was found
% thanks to a nice bradycardia of baby4 (cECG
% study\20120704_test4\20120704_182322)

%V13-4-2017
% made FS changeable in main_algorithm (was 250, has become 500) hence BUFFER_SIZE has
% changed from 1250 to 2500 to remain 5 seconds of data, 