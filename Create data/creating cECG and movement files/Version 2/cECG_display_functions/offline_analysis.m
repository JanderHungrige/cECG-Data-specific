TEXT_ONLY=1;
AVERAGE_GRAPHS=0;
timeSamp=[1:length(projectedECG)]./250./60+initOffset/60;
timeWin=[1:size(safeAmpli,2)]./8000*WINDOW_SIZE/60+initOffset/60;%+STARTtime; %in minutes %It has to be saved during the analysis IF the window sizes are dynamically changing. Here=asumed all same length.
position=Body*ones(size(timeWin));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializing text file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc
tic
cd(DataDirectory)
filename = sprintf('Results_%s.txt', DataDirectory(end-14:end)); % number serves only as an example
fid = fopen(filename,'a+t'); %w+t, open or create file (in text mo) for reading and writing; discard existing contents!

fprintf(fid,'%s\n','---------------------------------------------------------------------------------------------------------------------------');
fprintf(fid,'%s\t\t\t%s\n','Analysis date:', datestr(now,21)); % \t means "tab space"
fprintf(fid,'%s\t\t\t%s\n','Data set:', DataDirectory);
fprintf(fid,'%s\t\t%s\n','Sampling frequency:', num2str(FS));
fprintf(fid,'%s\t%s%s\n','Length processed signal:', num2str(sampleCounter/8000/60-initOffset/60), ' minutes');
fprintf(fid,'%s\t\t\t%s%s\n','Starting time:', num2str(initOffset/60), ' minutes after the beginning of the recording.');
if ALTERNATIVE_PEAK_DETECT
    fprintf(fid,'%s\n',['Configuration: WINDOW_SIZE=' num2str(WINDOW_SIZE)  ',KALMAN=' num2str(KALMAN) ',MIN3CHANNELS=' num2str(MIN3CHANNELS) ',CHAN_SELECT=' num2str(CHAN_SELECT) ',CHAN_2SELECT=' num2str(CHAN_2SELECT) ',FIXED_NUM_RR=' num2str(FIXED_NUM_RR) ',ANNOTATIONS=' num2str(ANNOTATIONS) 'RESP=' num2str(RESP) ',PRESSMAT=' num2str(PRESSMAT)]);
else
    fprintf(fid,'%s\n',['Configuration: WINDOW_SIZE=' num2str(WINDOW_SIZE)  ',KALMAN=' num2str(KALMAN) ',MIN3CHANNELS=' num2str(MIN3CHANNELS)  ',CHAN_SELECT=' num2str(CHAN_SELECT) ',FIXED_NUM_RR=' num2str(FIXED_NUM_RR) ',ANNOTATIONS=' num2str(ANNOTATIONS) 'RESP=' num2str(RESP) ',PRESSMAT=' num2str(PRESSMAT)]);

end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Motion-based segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% METHOD 1: RESP-based  (motionL)
% --------------------
if RESP
    THR1=5e4; %5e4 0.5e4 500
    THR2=1e6; %1e6 5e4 0.5e4
    fprintf(fid,'%s\t\t%g%s%g\n','Motion levels RESP Louis:', THR1, ' and ', THR2);
    safeMotion1kHz=abs(safeMotion1kHz);
    motionL=ones(size(safeMotion1kHz));
    motionL(safeMotion1kHz>=THR1)=2;
    motionL(safeMotion1kHz>=THR2)=3;
    if ~TEXT_ONLY
        f1=figure('Name', 'Motion from RESP Louis'); plot(timeWin,safeMotion1kHz);
        ax(6)=gca;
        plot_segments(f1,timeWin, motionL, min(safeMotion1kHz), abs(min(safeMotion1kHz))+max(safeMotion1kHz));
    end
    %matrixL=[timeWin',position',motionL',safeMotion1kHz',(safeMeanRR2')/250,(safeMeanRRREF')/250];
end

% METHOD 2: ampli-based  (motionA)
% --------------------

THR1=1e8; %5e8 for win size =4000 and winSize=5
THR2=1e9; %2e10 5e9
fprintf(fid,'%s\t\t%g%s%g\n','Motion levels Aline:', THR1, ' and ', THR2);
fprintf(fid,'%s\n','---------------------------------------------------------------------------------------------------------------------------');
safeAmpli2=sum(safeAmpli(:,:),1); %channelToUse online? No because no motion info!
%{
    % Method1: variance over last 2.5 seconds (online)
    motionIn=ones([1 size(safeAmpli,2)]);
    winSize=5;
    var_ampli=zeros(size(safeAmpli2));
    var_ampli(1:winSize-1)=var(safeAmpli2(1:winSize-1)); %give the same value to the winSize first samples (to make sure the signals data and data_filtered have the same length)
    for i=winSize:length(safeAmpli2)
        var_ampli(i)=var(safeAmpli2(i-winSize+1:i)); %compute the variance on a sliding windows. Here, the overlap is winSize-1.
        if var_ampli(i)>= THR1
            motionIn(i)=2;  % detect a huge peak when the variance is above a certain threshold
            if var_ampli(i)>= THR2
                motionIn(i)=3;
            end
        end
    end
%}
% Method2: variance over surrounding 2.5 seconds! (not causal)
motionA=ones([1 size(safeAmpli,2)]);
winSize=4;
var_ampli2=zeros(size(safeAmpli2));
var_ampli2(1:winSize/2)=var(safeAmpli2(1:winSize/2)); %give the same value to the winSize first samples (to make sure the signals data and data_filtered have the same length)
for i=winSize/2+1:length(safeAmpli2)-winSize/2
    var_ampli2(i)=var(safeAmpli2(i-winSize/2:i+winSize/2)); %compute the variance on a sliding windows. Here, the overlap is winSize-1.
    if var_ampli2(i)>= THR1
        motionA(i)=2;  % detect a huge peak when the variance is above a certain threshold
        if var_ampli2(i)>= THR2
            motionA(i)=3;
        end
    end
end

%{
    % Method 3: remove mean over past 5 s seconds and take sqrt
    motion3=ones([1 size(safeAmpli,2)]);
    winSize=10; %5 s
    sqrt_ampli=zeros(size(safeAmpli2));
    sqrt_ampli(1:winSize-1)=sqrt(abs(safeAmpli2(1:winSize-1)-mean(safeAmpli2(1:winSize-1)))); %give the same value to the winSize first samples (to make sure the signals data and data_filtered have the same length)
    for i=winSize:length(safeAmpli2)
        sqrt_ampli(i)=sqrt(abs(safeAmpli2(i)-mean(safeAmpli2(i-winSize+1:i)))); %compute the variance on a sliding windows. Here, the overlap is winSize-1.
        if sqrt_ampli(i)>= THR1
            motion3(i)=2;  % detect a huge peak when the variance is above a certain threshold
            if sqrt_ampli(i)>= THR2
                motion3(i)=3;
            end
        end
    end
%}

% %Method3: derivation of the safeAmpli
% motionA=ones([1 size(safeAmpli2,2)]);
% diff_ampli=diff(safeAmpli2,[],2);
% diff_ampli=[diff_ampli(:,1), diff_ampli]; %give the same value to the winSize first samples (to make sure the signals data and data_filtered have the same length)
% motionA(diff_ampli >= 1e8)=2;
% motionA(diff_ampli >= 1e9)=3;

% Plot the three motion levels
if ~TEXT_ONLY
    f1=figure('Name', 'Motion from 1kHz Aline'); plot(timeWin,var_ampli2);
    plot_segments(f1,timeWin, motionA, min(var_ampli2), abs(min(var_ampli2)+max(var_ampli2)));
    ax(11)=gca;
end

% METHOD 3: annotation-based  (motionANNO)
% --------------------------
if ANNOTATIONS
    annotations=xlsread([DataDirectory(1:end-30) 'Case Record Forms Clinical Study' DataDirectory(end-30:end-15) 'labels\' DataDirectory(end-14:end) '_labels.xlsx']);
    motionANNO=7*ones([1 length(timeWin)]);
    annotations(:,1:2)=annotations(:,1:2)-initOffset/60; %remove the offset (in minutes)
    annotations(:,1:2)=annotations(:,1:2).*60*8000/WINDOW_SIZE; % transform minutes into seconds then into samples of timeWin
    for i=1:size(annotations,1)
        motionANNO(max(0,annotations(i,1))+1:min(max(1,annotations(i,2)),length(timeWin)))=annotations(i,3); %annotations(1,1) should be = 0
    end
    if ~TEXT_ONLY
        f1=figure('Name', 'Motion from annotations'); plot(timeWin,zeros(size(timeWin)));
        plot_segments(f1,timeWin, motionANNO, -1, 2);
        ax(12)=gca;
    end
end

% METHOD 4: signal-based  (motionS)
% --------------------------
motionS=ones([1 size(projectedECG,2)]);
winSize=250; %250 = 1 second
var_ampli3=zeros(size(projectedECG));
if KALMAN
    var_ampli3(1:winSize/2)=var(cleanProjectedECG(1:winSize/2)); %give the same value to the winSize first samples (to make sure the signals data and data_filtered have the same length)
else
    var_ampli3(1:winSize/2)=var(projectedECG(1:winSize/2));
end
for i=winSize/2+1:length(projectedECG)-winSize/2
    var_ampli3(i)=var(projectedECG(i-winSize/2:i+winSize/2)); %compute the variance on a sliding windows. Here, the overlap is winSize-1.
    if var_ampli3(i)>= 6e6
        motionS(i)=2;  % detect a huge peak when the variance is above a certain threshold
        if var_ampli3(i)>= 12e6
            motionS(i)=3;
        end
    end
end
if ~TEXT_ONLY
    f1=figure('Name', 'Motion from signal itself'); plot(timeSamp, var_ampli3);
    plot_segments(f1,timeSamp, motionS, min(var_ampli3), abs(min(var_ampli3)+max(var_ampli3)));
    ax(13)=gca;
end


% METHOD 5: pressure-mat-based  (motionP)
% -----------------------------
if PRESSMAT
    THR1=2500;
    THR2=2e4;
    motionP=ones(size(safeMotionLevel));
    motionP(safeMotionLevel>=THR1)=2;
    motionP(safeMotionLevel>=THR2)=3;
    if ~TEXT_ONLY
        f1=figure('Name', 'Motion from PRESSURE MAT'); plot(timeWin,safeMotionLevel);
        ax(14)=gca;
        plot_segments(f1,timeWin, motionP, min(safeMotionLevel), abs(min(safeMotionLevel))+max(safeMotionLevel));
    end
end

% Plot all motion segmentation strategies
% ---------------------------------------
if ~TEXT_ONLY
    f1=figure('Name', 'Motion segmentation comparison');
    titlee='Segmentation (';
    if ANNOTATIONS
        plot_segments(f1, timeWin, motionANNO, 0, 1);
        titlee=[titlee 'Anno,'];
    end
    if RESP
        plot_segments(f1, timeWin, motionL, 2, 1);
        titlee=[titlee 'Louis,'];
    end
    if PRESSMAT
        plot_segments(f1, timeWin, motionP, -2, 1);
        titlee=[titlee 'Pressmat,'];
    end
    plot_segments(f1, timeWin, motionA, 1, 1);
    plot_segments(f1, timeSamp, motionS, -1, 1);
    ylim([-2 2])
    xlim([0 timeWin(end)])
    titlee=[titlee 'Aline,Signal)'];
    title(titlee)
    xlabel('Time [min]')
    ax(15)=gca;
end
% Decide which segmentation method to use for HR evaluation
% ---------------------------------------------------------
if ANNOTATIONS
    motionIn=motionANNO;fprintf(fid,'%s\n','Motion segmentation based on: annotations');
elseif RESP
    motionIn=motionL;fprintf(fid,'%s\n','Motion segmentation based on: 1 kHz (Louis)');
elseif PRESSMAT
    motionIn=motionP;fprintf(fid,'%s\n','Motion segmentation based on: pressure-mat');
else
    motionIn=motionA;fprintf(fid,'%s\n','Motion segmentation based on: 1 kHz (Aline)');
end

%matrixA=[timeWin',position',motionIn',(safeMeanRR')./FS,(safeMeanRR2')./FS,(safeMeanRRREF')./FS]; %var_ampli2 %RR intervals in seconds

% Upsample motion (from 2Hz (if WINDOW_SIZE=4000) to 250 Hz)
motionUp=repmat(motionIn,WINDOW_SIZE/32,1);
motionUp=motionUp(:)';
%NB: motionIn @8000/WINDOW_SIZE Hz = motionUp@250 Hz = motion @1 per peak position (unevenly sampled)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Averaged HR evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if AVERAGE_GRAPHS
    f1=figure('Name', 'Average HR'); plot(timeWin,60*FS./safeMeanRR);
    hold on; plot(timeWin,60*FS./safeMeanRRREF,'r')
    ax(9)= gca;
    plot_segments(f1,timeWin, motionIn, min(60*FS./safeMeanRRREF), abs(min(max(60*FS./safeMeanRRREF),HR(2))-min(60*FS./safeMeanRRREF)));
    f1=figure('Name', 'Average HR (Rel Ind based)'); plot(timeWin,60*FS./safeMeanRR2);
    hold on; plot(timeWin,60*FS./safeMeanRRREF,'r')
    plot_segments(f1,timeWin, motionIn, min(60*FS./safeMeanRRREF), abs(min(max(60*FS./safeMeanRRREF),HR(2))-min(60*FS./safeMeanRRREF)));
    ax(10)= gca;
end

% Compute the error for error/tolerance table
if ~FIXED_NUM_RR
    fprintf(fid,'%s%2.1f%s\n','Results when heart rate averaged over the last ', NUM_RR_AVERAGED  ,' seconds');
else
    fprintf(fid,'%s%2f%s\n','Results when heart rate averaged over the last ', NUM_RR_AVERAGED  ,' RR intervals');
end
fprintf(fid,'%s\n','---------------------------------------------------------------------------------------------------------------------------');
error0 = abs(60*FS./safeMeanRR-60*FS./safeMeanRRREF); % in bpm %1 value per processing windows: FS=8000/WINDOW_SIZE
fprintf(fid,'%s\t%s\t%s\t%s\n','Tolerance', 'MatchWithREF(%time)', 'NoMatch(#segments>10s)', 'MaxSegment(s)');
evaluate_matching(fid, error0, 1, WINDOW_SIZE); %tolerance of 1 beat per minute deviation with the ref
evaluate_matching(fid, error0, 3, WINDOW_SIZE); %3 bpm
evaluate_matching(fid, error0, 5, WINDOW_SIZE); %4 bpm
fprintf(fid, '\n');

fprintf(fid,'%s\n','Results when a reliability indicator is used to discard the unexpected RR intervals in the computation of the average:');
fprintf(fid, '\n');
errorRelInd = abs(60*FS./safeMeanRR2-60*FS./safeMeanRRREF); % in bpm
fprintf(fid,'%s\t%s\t%s\t%s\n','Tolerance', 'MatchWithREF(%time)', 'NoMatch(#segments>10s)', 'MaxSegment(s)');
saved_length=evaluate_matching(fid, errorRelInd, 1, WINDOW_SIZE); %tolerance of 1 beat per minute deviation with the ref
evaluate_matching(fid, errorRelInd, 3, WINDOW_SIZE); %3 bpm
evaluate_matching(fid, errorRelInd, 5, WINDOW_SIZE); %4 bpm
fprintf(fid, '\n');

if ~TEXT_ONLY
    % Histogram for 1 bpm tolerance WITH the reliability indicator
    figure('Name', 'Histogram (with rel_indic)');
    hist(saved_length); title('Segment lengths distribution for 1 bpm tolerance');
    xlabel('Segment length (s)')
    ylabel('Number of segments')
end

% Table with motion segmentation
fprintf(fid,'%s\n\n', 'Match with the reference (in % of time) for different motion levels and different deviation tolerances (and when reliability indicator used):');
fprintf(fid,'\t\t\t%s\t\t%s\t\t%s\t\t%s\n', 'Length(min)', '1 bpm', '3 bpm', '5 bpm');
% No motion
error_temp=errorRelInd;
fprintf(fid,'%s\t\t%3.2f\t\t\t%2i\t\t%2i\t\t%2i\n', 'Whole signal', length(motionIn)/8000*WINDOW_SIZE/60, round(sum(error_temp<=1)/length(error_temp)*100), round(sum(error_temp<=3)/length(error_temp)*100), round(sum(error_temp<=5)/length(error_temp)*100));
fprintf(fid, '\n');
error_temp=errorRelInd(motionIn==1);
fprintf(fid,'%s\t\t%3.2f\t\t\t%2i\t\t%2i\t\t%2i\n', 'No motion   ', sum(motionIn==1)/8000*WINDOW_SIZE/60, round(sum(error_temp<=1)/length(error_temp)*100), round(sum(error_temp<=3)/length(error_temp)*100), round(sum(error_temp<=5)/length(error_temp)*100));
error_temp=errorRelInd(motionIn==2);
fprintf(fid,'%s\t\t%3.2f\t\t\t%2i\t\t%2i\t\t%2i\n', 'Low motion  ', sum(motionIn==2)/8000*WINDOW_SIZE/60, round(sum(error_temp<=1)/length(error_temp)*100), round(sum(error_temp<=3)/length(error_temp)*100), round(sum(error_temp<=5)/length(error_temp)*100));
error_temp=errorRelInd(motionIn==3);
fprintf(fid,'%s\t\t%3.2f\t\t\t%2i\t\t%2i\t\t%2i\n', 'High motion ', sum(motionIn==3)/8000*WINDOW_SIZE/60, round(sum(error_temp<=1)/length(error_temp)*100), round(sum(error_temp<=3)/length(error_temp)*100), round(sum(error_temp<=5)/length(error_temp)*100));
error_temp=errorRelInd(motionIn==4|motionIn==5);
fprintf(fid,'%s\t\t%3.2f\t\t\t%2i\t\t%2i\t\t%2i\n', 'Intervention', sum(motionIn==4|motionIn==5)/8000*WINDOW_SIZE/60, round(sum(error_temp<=1)/length(error_temp)*100), round(sum(error_temp<=3)/length(error_temp)*100), round(sum(error_temp<=5)/length(error_temp)*100));
error_temp=errorRelInd(motionIn==6);
fprintf(fid,'%s\t\t%3.2f\t\t\t%2i\t\t%2i\t\t%2i\n', 'Ghosts      ', sum(motionIn==6)/8000*WINDOW_SIZE/60, round(sum(error_temp<=1)/length(error_temp)*100), round(sum(error_temp<=3)/length(error_temp)*100), round(sum(error_temp<=5)/length(error_temp)*100));
error_temp=errorRelInd(motionIn==7);
fprintf(fid,'%s\t\t%3.2f\t\t\t%2i\t\t%2i\t\t%2i\n', 'Ignored Data', sum(motionIn==7)/8000*WINDOW_SIZE/60, round(sum(error_temp<=1)/length(error_temp)*100), round(sum(error_temp<=3)/length(error_temp)*100), round(sum(error_temp<=5)/length(error_temp)*100));
fprintf(fid, '\n');

%NB: length of signal: here, the number of processing windows (WINDOW_SIZE)
%with the same motion label are summed to compute the length of the, say,
%no motion period. In the next table, each RR intervals are associated with
%a motion label (the one of the peak ending the interval) and the RR
%intervals are summed. Therefore there might be slight differencies.

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HR evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define artifacts on ref detections
intervals=diff(RRpositionsREF);
%unnaturals=find(intervals<=floor(60/HR(2)*FS)&intervals>=ceil(60/HR(1)*FS);
goodTime=sum(intervals(intervals>=floor(60/HR(2)*FS)&intervals<=ceil(60/HR(1)*FS)&intervals<1.5*mean(intervals)&intervals>0.5*mean(intervals)));
totalTime=sum(intervals);
disp(['Reasonable peak detection in the reference ECG during ' num2str(goodTime/totalTime*100) ' % of the time;'])
% hold on;
% temp=1.5*1/mean(intervals)*60*250;
% plot([1:RRpositionsREF(end)]./250,temp*ones(1,RRpositionsREF(end)/250))

%addpath('C:\Users\Aline\Documents\Research\Philips cECG\Dropbox\E-NEMO\Algorithm\rPeakDetection - 20111220michiel - PROTECTED\')
fprintf(fid,'%s\n','---------------------------------------------------------------------------------------------------------------------------');
fprintf(fid,'%s\n','Results on a peak to peak basis (intantaneous heart rate)');
fprintf(fid,'%s\n','---------------------------------------------------------------------------------------------------------------------------');
fprintf(fid,'%s\t\t%s\t%s\t%s\t%s\t%s\t\t\t%s\t%s\t%s\n','Signal    ', 'Length(min)', '%time',  '#seg>10s', 'Max(s)', 'Mean(s)', 'PPV', 'SEN', 'ERROR');
fprintf(fid, '\n');
[TPArray, FPArray, FNArray,results]=evaluate_detection_performance_per_motion_level(fid, RRpositions, RRpositionsREF, motionUp, FS, initOffset);%STARTtime);
%NB: motionIn @8000/WINDOW_SIZE Hz = motionUp@250 Hz = motion @1 per peak position (unevenly sampled)

if ALTERNATIVE_PEAK_DETECT
    fprintf(fid, '\n');
    fprintf(fid,'%s\n',['Results when ONLY the difference signal betwen TWO sensors is taken to detect the heart beats: sensor ' num2str(channelsToUseALT(1)) ' and sensor ' num2str(channelsToUseALT(2))]);
    fprintf(fid, '\n');
    fprintf(fid,'%s\t\t%s\t%s\t%s\t%s\t%s\t\t\t%s\t%s\t%s\n','Signal ALT',  'Length(min)', '%time',  '#seg>10s', 'Max(s)', 'Mean(s)', 'PPV', 'SEN', 'ERROR');
    fprintf(fid, '\n');
    [TPArrayALT, FPArrayALT, FNArrayALT,resultsALT]=evaluate_detection_performance_per_motion_level(fid, RRpositionsALT, RRpositionsREF, motionUp, FS, initOffset);
end

if ~TEXT_ONLY
    % histogram
    figure('Name', 'Histogram (iHR)');
    hist(results.wholeData(9:end)); title('Segment lengths distribution for iHR');
    xlabel('Segment length (s)')
    ylabel('Number of segments')
end
%}
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel selection analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~TEXT_ONLY
    %Signals:
    f1=figure('Name', 'Full graph'); plot(timeSamp, filteredECG(2,:)-filteredECG(4,:))
    ax(1)=gca;
    plot_segments(f1,timeWin, motionIn, -1e5, 1e5);
    hold on; plot(timeSamp, filteredECG(2,:)+20000,'g');
    hold on; plot(timeSamp, filteredECG(4,:)+15000,'c')
    %hold on; plot(timeSamp, filteredECG(2,:)-filteredECG(4,:),'b')
    hold on; plot(timeSamp, projectedECG-15000,'k')
    hold on; plot(timeSamp, filteredECG(9,:)./100-33000,'r')
    hold on; plot(RRpositions./250./60,projectedECG(RRpositions-round(initOffset*250))-15000,'om');
    fi=filteredECG(2,:)-filteredECG(4,:);
    hold on; plot(RRpositionsALT./250./60,fi(RRpositionsALT-round(initOffset*250)),'om');
    legend('Sen2 - Sen4', 'Sen2', 'Sen4', 'projected VCG', 'REF')
    %hold on; plot(RRpositionsREF./250./60,filteredECG(9,RRpositionsREF)./200-7000,'om');
    title('ECG signals')
    
    %hold on; plot([1:length(Rel_Indicator)]./8000.*numberNewSamples./60+STARTtime,Rel_Indicator(1,:).*(-50000)-43000,'.g')
    hold on; plot(timeWin,Rel_Indicator(1,:).*(-50000)-43000,'.g')
    
    hold on; plot(TPArray./250./60,ones(1,length(TPArray))-48000,'.r')
    hold on; plot(FPArray./250./60,ones(1,length(FPArray))-48000,'*r')
    hold on; plot(FNArray./250./60,ones(1,length(FNArray))-48000,'or')
    
    %Heart Rate:
    figure;
    if ALTERNATIVE_PEAK_DETECT
        stairs(RRpositionsALT(1:end-1)./250./60,250*60./diff(RRpositionsALT),'b');
    end
    hold on ; stairs(RRpositions(1:end-1)./250./60,250*60./diff(RRpositions),'k');
    hold on ; stairs(RRpositionsREF(1:end-1)./250./60,250*60./diff(RRpositionsREF),'r');
    ylim([90 220])
    ax(2)=gca;
    legend('iHR from Sen2 - Sen4', 'iHR from projected VCG', 'iHR from REF')
    title('Instantaneous heart rate')
    
    % %Channel selection:
    % figure; plot([1:length(safeSelect)]./8000.*numberNewSamples./60+STARTtime,safeSelect(1:8,:)+repmat([0.05 0.06 0.07 0.01 0.02 0.03 0.04 0.08]',1,size(safeSelect,2)),'.', 'MarkerSize',4)
    % ylim([1 1.1])
    % if ~ADULT
    %     legend('5','6','7','1','2','3','4','8')
    % else
    % end
    % ax(3)=gca;
    % title('Selected channels')
    figure; %plot([1:length(safeAmpli)]./8000.*numberNewSamples./60+STARTtime,safeAmpli(1:8,:))
    plot(timeWin,safeAmpli(1:8,:))
    % legend('1','2','3','4','5','6','7','8')
    % ax(4)=gca;
    % hold on ; plot([1:length(safeAmpli)]./8000.*numberNewSamples./60+STARTtime,1*median(safeAmpli), 'y')
    
    %Video-based motion level estimation
    if CAMERA
        figure; plot(timeWin,safeMotionLevel)
        ax(5)=gca;
    end
    
    
    
    % if RESP
    % %     %Respiration from 1kHz
    % %     %figure; plot([1:length(safeResp)]./250./60+initOffset/60,safeResp);
    % %     figure; plot(timeSamp,safeResp);
    % %     ax(7)=gca;
    % %     %Reference respiration
    % %     figure; plot(timeSamp,safeRespREF);
    % %     ax(8)=gca;
    % end
    
    linkaxes(ax,'x');
end
toc
%}

%%

%Reliability indicator

% sortedCombo(2,Rel_Indicator(min(ceil(sortedCombo(2,:)/(numberNewSamples/32)),numel(Rel_Indicator)))==sortedCombo(1,:));
% sortedCombo
% hold on; plot(sortedCombo(2,:)./250,sortedCombo(1,:).*2+48,'.r')
% time(1)=Rel_Indicator(2,1);
% for i=2:length(Rel_Indicator)
%     time(i)=time(i-1)+Rel_Indicator(2,i);
% end
% hold on; plot(time./250,Rel_Indicator(1,:).*2+43,'.g')
%timeStartECG
%hold on; plot([1:length(Rel_Indicator)].*62./250,Rel_Indicator.*2+43,'.b')
%figure; plot(Rel_Indicator(2,:)); %is it a constant?
%figure;pie([availableHRTime, totalTime-availableHRTime],[{'Correct iHR'},{'No iHR'}]);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
if ~GRAPH
    %toc
    %maxTIme
    %profile viewer
    %     meany=mean(profileARRAY,1);
    %     maxy=max(profileARRAY,[],1);
    %     miny=min(profileARRAY,[],1);
    %     figure; bar([meany(1:8)', maxy(1:8)']);
    %     title('Execution time for processing windows of 1 second')
    %     legend('Mean execution time', 'Maximum execution time')
    %     ylabel('Time (s)')
    %     set(gca,'XTickLabel',{'Chan.','Down.','Bipol.','Band.','VCG','Peaks','Kalman','TOTAL'},'XTick',[1 2 3 4 5 6 7 8]);
end
if ~isempty(RRpositions)
    figure;
    if numel(RRpositionsREF)>1
        stairs(RRpositionsREF(1:end-1)./FS_REF,60*FS_REF./diff(RRpositionsREF),'r'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
    end
    hold on;
    stairs(RRpositions(1:end-1)./FS,60*FS./diff(RRpositions),'b'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
    
    %    stairs(RRpositionsALT(1:end-1)./FS,60*FS./diff(RRpositionsALT),'g'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
    %legend('with channel selection', 'with channel selection and reliability test', 'Reference')
    
    
    
    %Run Michiel's code on the dataset for comparison purpose
    %     addpath([cd '\rPeakDetection - 20111220michiel - PROTECTED\'])
    %     [out] = streamingpeakdetection(projectedECG, FS, HR, false, 16); %HRLimits, PLOT, Fc, blockSize)
    %     peakPositionArray1 = out.peakPositionArray;
    %     hold on;
    %     stairs(peakPositionArray1(1:end-1)./FS,60*FS./diff(peakPositionArray1),'c'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
    
    %figure; stairs([1:numel(Rel_Indicator)].*numberNewSamples,Rel_Indicator,'k')
    
    if FS_REF~=8000
        %        [out] = streamingpeakdetection(referenceECG, FS_REF); %HRLimits, PLOT, Fc, blockSize)
        %        peakPositionArray = out.peakPositionArray;
        %        hold on;
        %stairs(peakPositionArray(1:end-1)./FS_REF,60*FS_REF./diff(peakPositionArray),'m'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
    end
    
    
    % perform peak detection on tempData
    
    %addpath('C:\Users\Aline\Documents\Research\Apnea\MATLAB')
    %line_list=swab(tempDataECG, 'fs', 250,'graph',false);
    %allPeakLocation=line_list(2:2:end-1,2);
    allPeakLocation=RRpositions';
    allPeakLocationALT=RRpositionsALT';
    
    %     % filter
    %     RRintervals=diff(allPeakLocation);
    %     indexes=find(RRintervals<=60);%250 bpm
    %     % fuse two consecutive small RRintervals
    %     indexes=indexes(indexes(1:end-1)==indexes(2:end)-1);
    %     RRintervals(indexes)=RRintervals(indexes)+RRintervals(indexes+1);
    %     RRintervals(indexes+1)=[];
    %     allPeakLocation=[allPeakLocation(1); cumsum(RRintervals)+allPeakLocation(1)];
    %     %hold on;
    %     %stairs(allPeakLocation(1:end-1)./FS,60*FS./diff(allPeakLocation),'k')
    
    % filter
    %     peakToRemove=[];
    %     for i=1:numel(allPeakLocation)
    %         index=find(allPeakLocationALT>=allPeakLocation(i)-25 & allPeakLocationALT<=allPeakLocation(i)+25); %100ms
    %         if isempty(index)
    %             peakToRemove=[peakToRemove i];
    %         end
    %     end
    %     allPeakLocation(peakToRemove)=[];
    RRintervals=diff(allPeakLocation);
    
    %figure;
    RRintervalsMean=[];
    for i=16:numel(RRintervals)
        RRintervalsMean(i-15)=mean(RRintervals(i-15:i));
    end
    %stairs(allPeakLocation(16:end-1)./FS,60*FS./ RRintervalsMean,'b')
    
    
    RRintervals2=diff(allPeakLocationALT);
    RRintervalsMean=[];
    for i=16:numel(RRintervals2)
        RRintervalsMean(i-15)=mean(RRintervals2(i-15:i));
    end
    %    hold on; stairs(allPeakLocationALT(16:end-1)./FS,60*FS./ RRintervalsMean,'g')
    
    %     RRintervalsMean=[];
    %     RRintervals=diff(peakPositionArray);
    %     for i=16:numel(RRintervals)
    %         RRintervalsMean(i-15)=mean(RRintervals(i-15:i));
    %     end
    %     hold on; stairs(peakPositionArray(16:end-1)./FS,60*FS./RRintervalsMean,'r')
    
    % figure; plot([1:numel(projectedECG)]./FS, projectedECG)
end


%plotdata(8,[1:size(filteredECG, 2)]/FS,filteredECG,'NOISY SIGNAL (time)');

%}

%%%%%%%%%%%%%%%
% GHOST BUSTER
%%%%%%%%%%%%%%%

%{

% No try
figure;ax(1)=subplot(3,1,1);plot(projectedECG);
title('NO TRY: projected ECG')
ax(2)=subplot(3,1,2);plot(filteredECG(4,:)-filteredECG(1,:));
title('sens 4 - sens 1')
ax(3)=subplot(3,1,3);plot(filteredECG(5,:)-filteredECG(6,:));
title('sens 5 - sens 6')
linkaxes(ax,'x')

figure; ax(1)=subplot(3,1,1);
stairs(RRpositionsREF(1:end-1)./FS_REF,60*FS_REF./diff(RRpositionsREF),'r'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
hold on;
peakLoca=offlinePeakDetect(projectedECG,filterCoefficients);
stairs(peakLoca(1:end-1)./FS,60*FS./diff(peakLoca),'b'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
title('NO TRY')
ax(2)=subplot(3,1,2);
stairs(RRpositionsREF(1:end-1)./FS_REF,60*FS_REF./diff(RRpositionsREF),'r'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
hold on;
peakLoca=offlinePeakDetect(filteredECG(4,:)-filteredECG(1,:),filterCoefficients);
stairs(peakLoca(1:end-1)./FS,60*FS./diff(peakLoca),'b'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
title('sens 4 - sens 1')
ax(3)=subplot(3,1,3);
stairs(RRpositionsREF(1:end-1)./FS_REF,60*FS_REF./diff(RRpositionsREF),'r'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
hold on;
peakLoca=offlinePeakDetect(filteredECG(5,:)-filteredECG(6,:),filterCoefficients);
stairs(peakLoca(1:end-1)./FS,60*FS./diff(peakLoca),'b'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
title('sens 5 - sens 6')
linkaxes(ax,'x')

%input('now what?')

%%%%%%%%%%%%%%%%%%%%%%%%% 5 and 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x1=filteredECG(5,1:14730);
% x2=filteredECG(6,1:14730);
%
% L=size(x1,2);
% NFFT = 2^nextpow2(L); % Next power of 2 from length of signal
% f = FS/2*linspace(0,1,NFFT/2+1);
% F1=fft(x1,NFFT)/L;
% F2=fft(x2,NFFT)/L;
% figure; hha(1)=subplot(2,1,1);plot(f,abs(F1(1:NFFT/2+1))./abs(F2(1:NFFT/2+1)))
% title('Ratio of amplitudes between sens5 and sens6 on ECG')
% hha(2)=subplot(2,1,2);plot(f,angle(F1(1:NFFT/2+1))-angle(F2(1:NFFT/2+1)))
% title('Difference of phases between sens5 and sens6 on ECG')
% linkaxes(hha,'x')
% xlim([3 50])
%ylim([-0.5 0.5]);
%%

x1=filteredECG(6,:);
x2=filteredECG(5,:);

L=size(x1,2);
NFFT = 2^nextpow2(L); % Next power of 2 from length of signal
f = FS/2*linspace(0,1,NFFT/2+1);
F1=fft(x1,NFFT)/L;
F2=fft(x2,NFFT)/L;
figure; hha(1)=subplot(2,1,1);plot(f,abs(F1(1:NFFT/2+1))./abs(F2(1:NFFT/2+1)))
hold on
title('Ratio of amplitudes between sens6 and sens5')
hha(2)=subplot(2,1,2);plot(f,angle(F1(1:NFFT/2+1))-angle(F2(1:NFFT/2+1)))
hold on
title('Difference of phases between sens6 and sens5')
linkaxes(hha,'x')
xlim([3 50])
%ylim([-0.5 0.5]);

fc=0.5;
RC=1/(2*pi*fc);
alphaa=RC/(RC+1/250);
filtx1=0;
for i=2:size(x1,2)
    filtx1(i)=alphaa*(filtx1(i-1)+x1(i)-x1(i-1));
end

x1=filtx1;
x2=x2;

L=size(x1,2);
NFFT = 2^nextpow2(L); % Next power of 2 from length of signal
f = FS/2*linspace(0,1,NFFT/2+1);
F1=fft(x1,NFFT)/L;
F2=fft(x2,NFFT)/L;
plot(hha(1),f,abs(F1(1:NFFT/2+1))./abs(F2(1:NFFT/2+1)),'g')
title('Ratio of amplitudes between filt sens6 and sens5')
plot(hha(2),f,angle(F1(1:NFFT/2+1))-angle(F2(1:NFFT/2+1)),'g')
title('Difference of phases between filt sens6 and sens5')
linkaxes(hha,'x')
xlabel('Frequency (Hz)')
xlim([3 50])

figure;ax(1)=subplot(3,1,1);plot(projectedECG);
title('projected ECG')
ax(2)=subplot(3,1,2);plot(filteredECG(6,:)-filteredECG(5,:));
title('sens 6 - sens 5')
ax(3)=subplot(3,1,3);plot(filtx1-filteredECG(5,:),'g');
title('filtsens 6 - sens 5')
linkaxes(ax,'x')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SEN1 SEN4 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
x1=filteredECG(4,:);
x2=filteredECG(1,:);

L=size(x1,2);
NFFT = 2^nextpow2(L); % Next power of 2 from length of signal
f = FS/2*linspace(0,1,NFFT/2+1);
F1=fft(x1,NFFT)/L;
F2=fft(x2,NFFT)/L;
figure; hha(1)=subplot(2,1,1);plot(f,abs(F1(1:NFFT/2+1))./abs(F2(1:NFFT/2+1)))
hold on
title('Ratio of amplitudes between sens4 and sens1 on GHOST+ECG')
hha(2)=subplot(2,1,2);plot(f,angle(F1(1:NFFT/2+1))-angle(F2(1:NFFT/2+1)))
hold on
title('Difference of phases between sens4 and sens1 on GHOST+ECG')
linkaxes(hha,'x')
xlim([3 50])

%Filter x1 to compensate for the unbalance
% 1. manualy
fc=0.81;
RC=1/(2*pi*fc);
alphaa=RC/(RC+1/250)
filtx1=0;
for i=2:size(x1,2)
    filtx1(i)=alphaa*(filtx1(i-1)+x1(i)-x1(i-1));
end
%fircls1(n,wo,dp,ds)

x1=filtx1;
x2=x2;

L=size(x1,2);
NFFT = 2^nextpow2(L); % Next power of 2 from length of signal
f = FS/2*linspace(0,1,NFFT/2+1);
F1=fft(x1,NFFT)/L;
F2=fft(x2,NFFT)/L;
plot(hha(1),f,abs(F1(1:NFFT/2+1))./abs(F2(1:NFFT/2+1)),'g')
title('Ratio of amplitudes between filt sens4 and sens1')
plot(hha(2),f,angle(F1(1:NFFT/2+1))-angle(F2(1:NFFT/2+1)),'g')
title('Difference of phases between filt sens4 and sens1')
linkaxes(hha,'x')
xlabel('Frequency (Hz)')
xlim([3 50])

figure;ax(1)=subplot(3,1,1);plot(projectedECG(:));
title('NO TRY: projected ECG')
ax(2)=subplot(3,1,2);plot(filteredECG(4,:)-filteredECG(1,:));
title('sens 4 - sens 1')
ax(3)=subplot(3,1,3);plot(filtx1-filteredECG(1,:));
title('filtsens 4 - sens 1')
linkaxes(ax,'x')
%%
% f = FS/2*linspace(0,1,NFFT/2+1);
% %spectra = zeros(8,numel(f));
% for i = 1:size(data,1)
%     Y = fft(data(i,:),NFFT)/size(data,2);
%     %spectra(i,:) = 2*abs(Y(1:NFFT/2+1));
% end

%
% % Try1: adaptive filter
% [filterOutput, coefficients, weightedMean]=adaptiveFilter(filteredECG(1:8,:),mean(filteredECG(1:8,:)),250,2,'averageUpdate'); %coefficients of channel 1
% figure;ax(1)=subplot(3,1,1);plot(filteredECG(6,:)-filterOutput(6,:));
% title('TRY 1: adaptive filter with updated mean (sens 6 - sensFilt 6')
% ax(2)=subplot(3,1,2);plot(filterOutput(4,:)-filterOutput(1,:));
% title('TRY 1: sensFilt 4 - sensFilt 1')
% ax(3)=subplot(3,1,3);plot(filterOutput(5,:)-filterOutput(6,:));
% title('TRY1: sensFilt 5 - sensFilt 6')
% linkaxes(ax,'x')
%
% figure; ax(1)=subplot(3,1,1);
% stairs(RRpositionsREF(1:end-1)./FS_REF,60*FS_REF./diff(RRpositionsREF),'r'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
% hold on;
% peakLoca=offlinePeakDetect(filteredECG(6,:)-filterOutput(6,:),filterCoefficients);
% stairs(peakLoca(1:end-1)./FS,60*FS./diff(peakLoca),'b'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
% title('TRY 1: adaptive filter with updated mean (sens 6 - sensFilt 6')
% ax(2)=subplot(3,1,2);
% stairs(RRpositionsREF(1:end-1)./FS_REF,60*FS_REF./diff(RRpositionsREF),'r'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
% hold on;
% peakLoca=offlinePeakDetect(filterOutput(4,:)-filterOutput(1,:),filterCoefficients);
% stairs(peakLoca(1:end-1)./FS,60*FS./diff(peakLoca),'b'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
% title('TRY 1: sensFilt 4 - sensFilt 1')
% ax(3)=subplot(3,1,3);
% stairs(RRpositionsREF(1:end-1)./FS_REF,60*FS_REF./diff(RRpositionsREF),'r'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
% hold on;
% peakLoca=offlinePeakDetect(filterOutput(5,:)-filterOutput(6,:),filterCoefficients);
% stairs(peakLoca(1:end-1)./FS,60*FS./diff(peakLoca),'b'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
% title('TRY 1: sensFilt 5 - sensFilt 6')
% linkaxes(ax,'x')

% Try2: lowpass filter then adaptative
load('lowlowpass3_Fs250'); delay=116;
%load('lowpass_250_5_12.mat') %4 12 lowpass 250, 90 Num
%coefficients=Num; delay=45;
%artifact=filter(coefficients,1,mean(filteredECG(1:8,:)),[],2);
artifact=filter(coefficients,1,projectedECG,[],2);
[filterOutput, junk, weightedMean]=adaptiveFilter(artifact(delay:end),projectedECG(1:end-delay+1),250,1,'simple');

figure;ax(1)=subplot(3,1,1);plot(projectedECG(1:end-delay+1));
title('projected ECG')
ax(2)=subplot(3,1,2);plot(projectedECG(1:end-delay+1)-artifact(delay:end));
title('TRY 2: projected ECG - artifact (LPF)')
ax(3)=subplot(3,1,3);plot(projectedECG(1:end-delay+1)-filterOutput);
title('TRY 2: projected ECG - artifactFilt (LPF)')
linkaxes(ax,'x')

figure; ax(1)=subplot(3,1,1);
stairs(RRpositionsREF(1:end-1)./FS_REF,60*FS_REF./diff(RRpositionsREF),'r'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
hold on;
peakLoca=offlinePeakDetect(projectedECG(1:end-delay+1),filterCoefficients);
stairs(peakLoca(1:end-1)./FS,60*FS./diff(peakLoca),'b'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
title('projected ECG')
ax(2)=subplot(3,1,2);
stairs(RRpositionsREF(1:end-1)./FS_REF,60*FS_REF./diff(RRpositionsREF),'r'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
hold on;
peakLoca=offlinePeakDetect(projectedECG(1:end-delay+1)-artifact(delay:end),filterCoefficients);
stairs(peakLoca(1:end-1)./FS,60*FS./diff(peakLoca),'b'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
title('TRY 2: projected ECG - artifact (LPF)')
ax(3)=subplot(3,1,3);
stairs(RRpositionsREF(1:end-1)./FS_REF,60*FS_REF./diff(RRpositionsREF),'r'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
hold on;
peakLoca=offlinePeakDetect(projectedECG(1:end-delay+1)-filterOutput,filterCoefficients);
stairs(peakLoca(1:end-1)./FS,60*FS./diff(peakLoca),'b'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
title('TRY 2: projected ECG - artifactFilt (LPF)')
linkaxes(ax,'x')
a=1;


% Try3: median filter then adaptative
artifact=medianFilter(projectedECG,125); %0.5s?
[filterOutput, coefficients, weightedMean]=adaptiveFilter(artifact,projectedECG,250,2,'simple');

figure;ax(1)=subplot(3,1,1);plot(projectedECG(1:end));
title('projected ECG')
ax(2)=subplot(3,1,2);plot(projectedECG-artifact);
title('TRY 3: projected ECG - artifact (median)')
ax(3)=subplot(3,1,3);plot(projectedECG-filterOutput);
title('TRY 3: projected ECG - artifactFilt (median)')
linkaxes(ax,'x')


figure; ax(1)=subplot(3,1,1);
stairs(RRpositionsREF(1:end-1)./FS_REF,60*FS_REF./diff(RRpositionsREF),'r'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
hold on;
peakLoca=offlinePeakDetect(projectedECG,filterCoefficients);
stairs(peakLoca(1:end-1)./FS,60*FS./diff(peakLoca),'b'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
title('projected ECG')
ax(2)=subplot(3,1,2);
stairs(RRpositionsREF(1:end-1)./FS_REF,60*FS_REF./diff(RRpositionsREF),'r'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
hold on;
peakLoca=offlinePeakDetect(projectedECG-artifact,filterCoefficients);
stairs(peakLoca(1:end-1)./FS,60*FS./diff(peakLoca),'b'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
title('TRY 3: projected ECG - artifact (median)')
ax(3)=subplot(3,1,3);
stairs(RRpositionsREF(1:end-1)./FS_REF,60*FS_REF./diff(RRpositionsREF),'r'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
hold on;
peakLoca=offlinePeakDetect(projectedECG-filterOutput,filterCoefficients);
stairs(peakLoca(1:end-1)./FS,60*FS./diff(peakLoca),'b'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
title('TRY 3: projected ECG - artifactFilt (median)')
linkaxes(ax,'x')
a=1;



%Try4
% Try2: lowpass filter then adaptative
load('lowlowpass3_Fs250'); delay=116;
%load('lowpass_250_5_12.mat') %4 12 lowpass 250, 90 Num
%coefficients=Num; delay=45;
artifact=filter(coefficients,1,mean(filteredECG(1:8,:)),[],2);
[filterOutput, junk, weightedMean]=adaptiveFilter(artifact(delay:end),projectedECG(1:end-delay+1),250,1,'simple');

figure;ax(1)=subplot(3,1,1);plot(projectedECG(1:end-delay+1));
title('projected ECG')
ax(2)=subplot(3,1,2);plot(projectedECG(1:end-delay+1)-artifact(delay:end));
title('TRY 2: projected ECG - artifact (LPF)')
ax(3)=subplot(3,1,3);plot(projectedECG(1:end-delay+1)-filterOutput);
title('TRY 2: projected ECG - artifactFilt (LPF)')
linkaxes(ax,'x')

figure; ax(1)=subplot(3,1,1);
stairs(RRpositionsREF(1:end-1)./FS_REF,60*FS_REF./diff(RRpositionsREF),'r'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
hold on;
peakLoca=offlinePeakDetect(projectedECG(1:end-delay+1),filterCoefficients);
stairs(peakLoca(1:end-1)./FS,60*FS./diff(peakLoca),'b'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
title('projected ECG')
ax(2)=subplot(3,1,2);
stairs(RRpositionsREF(1:end-1)./FS_REF,60*FS_REF./diff(RRpositionsREF),'r'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
hold on;
peakLoca=offlinePeakDetect(projectedECG(1:end-delay+1)-artifact(delay:end),filterCoefficients);
stairs(peakLoca(1:end-1)./FS,60*FS./diff(peakLoca),'b'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
title('TRY 4: projected ECG - artifact (LPFmean)')
ax(3)=subplot(3,1,3);
stairs(RRpositionsREF(1:end-1)./FS_REF,60*FS_REF./diff(RRpositionsREF),'r'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
hold on;
peakLoca=offlinePeakDetect(projectedECG(1:end-delay+1)-filterOutput,filterCoefficients);
stairs(peakLoca(1:end-1)./FS,60*FS./diff(peakLoca),'b'); title('Instantaneous heart rate'); xlabel('Time (s)'); ylabel('HR (bpm)')
title('TRY 4: projected ECG - artifactFilt (LPFmean)')
linkaxes(ax,'x')
a=1;
% compute here your oracle from 5-6, 4-1, etc

%}

fclose(fid);
cd(CurrentDirectory)


