function [cECG1, cECG2, HR, motionA, motionL, Nr_channels_used, raw_cECG, safeResp] = cECG_extraction(raw_signals)
%
% [cECG1, cECG2, HR, motionA, motionL, Nr_channels_used, raw_cECG] = cECG_extraction(raw_signals)
%
% INPUT
% -----
%
% - raw_signals: 8xT matrix of each of the 8 unipolar ECG leads from the capacitive
% sensor array
%
%
% OUTPUT
% ------
%
% - cECG1:  1 bipolar ECG lead (EInthoven lead II) resulting from the selection and fusion of the 8 raw
% unipolar signals: length T/8000*500, Fs = 500 Hz
%
% - cECG2: 1 bipolar ECG lead  resulting from the simple difference between
% two channels (2 and 4): length T/8000*500, Fs = 500 Hz
%
% - HR: mean heart rate: length T/8000, Fs = 1 Hz as bpm
%
% - motionA: contact signals1KHz SIgnal injected. a measure of the amount of motion: length = T/8000, Fs = 1 Hz
%
% - motionL: another measure of the amount of motion: length = T/8000, Fs = 1 Hz
%
% - Nr_channels_used: vector that contains the number of capacitive channels (1-8) used for the extraction of cECG1. length = T/8000, Fs = 1 Hz
%
% - raw_cECG: the 8 capacitive channels (8 unipolar ECGs) with only band-pass
% filtering and downsampling: length T/8000*500, Fs = 500 Hz
%
%
%
%
% 26-06-2017
% made for Jan by aline.serteyn@philips.com

addpath([cd '\cECG_processing_functions'])
addpath([cd '\cECG_display_functions'])

% Fixed values:
% ------------
ONLINE_DEMO=0;
DO_INIT=1;
KALMAN_RESET=1;
PEAK_RESET=1;
GRAPH=0;
CAMERA=0;
PRESSMAT=0; %read the pressure mat files for motion level estimation
ANNOTATIONS=0; %read the annotations (.xlsx files) for motion level estimation
ALTERNATIVE_PEAK_DETECT=0; %if ALTERNATIVE_PEAK_DETECT=1 --> detect peaks on an alternative signal (e.g. diff of 2 fixed sensors) in parallel with the normal processing (used for offline analysis)
CHAN_2SELECT=0; %WHEN ALTERNATIVE_PEAK_DETECT=1 --> CHAN_SELECT=0=always #2-#4; CHAN_SELECT=1=always 2 best coupled channels; CHAN_SELECT=2=always 2 best non-adjacents
CHAN_SELECT=0; %CHAN_SELECT=0=old Aline selection strategy, CHAN_SELCT=1=new Louis selection strategy (TO BE IMPLEMENTED in main_algorithm.m (line303))
Manual=1; %body position determination (always 1 until the automated one is implemented)
SAVE_TO_DIR = 0; %save the cap ecg and HR to a given directory. This is done per dataset, hence per recording session (per data folder).
    
% May be changed:
% --------------
VIDEO=0; %1 may work for short demo time (reasonable memory usage)
RESP=1;
KALMAN=0; %apply Kalman smoother to the data (1=yes, 0=no)
VERBOSE=0;
SAVE=1; %1 may work for short demo time (reasonable memory usage)
ADULT=0; %0 for neonate, 1 for adults
Body = 1; %1 for chest, 0 for back
NOTCH24=1; %if =1: applies a 24Hz notch filter to the data
FIXED_NUM_RR=0; %choice of averaging strategy for the average heart rate (meanRR): 1=fixed number of RRintervals used for the computation of the meanRR, 0=last 10seconds (at least) of data used for the computation of the meanRR
MIN3CHANNELS=1; %force the selection of minimum 3 covered (not necessarly nicely coupled) channels in order to compute the 3 Einthoven leads

WIN_SIZE = 8000;
WINDOW_SIZE=8000;

for i = WIN_SIZE:WIN_SIZE:length(raw_signals)
    ECG1 = raw_signals(8,i-WIN_SIZE+1:i); %no reference ECG available sojust give a random capacitive channel
    SEN1 = raw_signals(1,i-WIN_SIZE+1:i);
    SEN2 = raw_signals(2,i-WIN_SIZE+1:i);
    SEN3 = raw_signals(3,i-WIN_SIZE+1:i);
    SEN4 = raw_signals(4,i-WIN_SIZE+1:i);
    SEN5 = raw_signals(5,i-WIN_SIZE+1:i);
    SEN6 = raw_signals(6,i-WIN_SIZE+1:i);
    SEN7 = raw_signals(7,i-WIN_SIZE+1:i);
    SEN8 = raw_signals(8,i-WIN_SIZE+1:i);
    if RESP
        referenceRespi=ECG1;
    end
    % !!!!!!!!!!!!!!!!!!!!!!!!!!
    main_algorithm
    % !!!!!!!!!!!!!!!!!!!!!!!!!!
end

% OUTPUTS
% -------
HR = 60*FS./safeMeanRR;
cECG1 = projectedECG;
cECG2 = filteredECG(2,:) - filteredECG(4,:);
Nr_channels_used = sum(safeSelect,1);
raw_cECG = filteredECG;

% Motion 
% ------
% Method1: variance over last 2.5 seconds (hence DELAY of 1.2 seconds!!)
% THR1=1e8; %5e8 for win size =4000 and winSize=5
% THR2=1e9; %2e10 5e9
safeAmpli2=sum(safeAmpli(:,:),1); 
winSize=5;
motionA=zeros(size(safeAmpli2));
motionA(1:winSize-1)=var(safeAmpli2(1:winSize-1)); %give the same value to the winSize first samples (to make sure the signals data and data_filtered have the same length)
for i=winSize:length(safeAmpli2)
    motionA(i)=var(safeAmpli2(i-winSize+1:i)); %compute the variance on a sliding windows. Here, the overlap is winSize-1.
end
% Motion: Method2
% THR1=5e4; %5e4 0.5e4 500
% THR2=1e6; %1e6 5e4 0.5e4
if RESP
    motionL = abs(safeMotion1kHz);
else
    motionL = NaN;
end


