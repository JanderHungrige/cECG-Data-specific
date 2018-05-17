%Create RR signal from DAQ and/or Intellivue signal
function RR=R_peak_detection_from_annotated_ECG(ECGsignal,FS_ecg,pat,padding,plotting)
    clc 
    clear
    tic

    Ralphsfactor={1;-1;-1;1;1;-1;-1;-1;1;1;-1;1;-1;1;-1;-1;-1;-1};%Determine if the ECG signal should be turned -1 or not 1. 

    Matlabfolder=('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific');
    addpath('C:\Users\310122653\Documents\PhD\Matlab\R peak detection and HRV\RAlps Rpeak detection')
    path='E:\';


    %%creating RR signal           
    t_DAQ=linspace(0,length(ECGsignal)/FS_ecg,length(ECGsignal))';
    RR=nan(1,length(ECG_bin_synched));
    [DAQ_RR_idx, ~, ~, ~, ~, DAQ_RR, ~] = ecg_find_rpeaks(t_DAQ,Ralphsfactor{Neonate,1}*ECGsignal, FS_ecg, 250,plotting,0); %, , , maxrate,plotting,saving   -1* because Ralph optimized for a step s slope, we also have steep Q slope. inverting fixes that probel 
    figure;title('RR')
    RR(DAQ_RR_idx)=DAQ_RR;

    disp('Finished')
    beep
end