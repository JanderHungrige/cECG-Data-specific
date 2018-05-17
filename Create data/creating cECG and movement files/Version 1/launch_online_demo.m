%LAUNCH_ONLINE_DEMO.M

% This is the initialization file for the online demo
% Simply run this in Matlab before starting the demo and call
% "main_algorithm" in the LabView/Matlab interface

addpath([cd '\cECG_processing_functions'])
addpath([cd '\cECG_display_functions'])

% Fixed values:
% ------------
ONLINE_DEMO=1;
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
    

% May be changed:
% --------------
VIDEO=0; %1 may work for short demo time (reasonable memory usage) 
RESP=1; 
KALMAN=1; %apply Kalman smoother to the data (1=yes, 0=no)
VERBOSE=0; 
SAVE=0; %1 may work for short demo time (reasonable memory usage) 
ADULT=0; %0 for neonate, 1 for adults
FIXED_NUM_RR=0; %choice of averaging strategy for the average heart rate (meanRR): 1=fixed number of RRintervals used for the computation of the meanRR, 0=last 10seconds (at least) of data used for the computation of the meanRR
MIN3CHANNELS=0; %force the selection of minimum 3 covered (not necessarly nicely coupled) channels in order to compute the 3 Einthoven leads
