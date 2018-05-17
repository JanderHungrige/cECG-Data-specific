% The EDR signal was all created with an RR detecion that was not optimal
% Therefore we create another one here with the correct ersion.


%Steps
% Load the synched ECG singal
% Create EDR signal with 1 or -1 as factor for the RR peak detection
% Save a s the new synched EDR
% 
% ame with Intellivuejust use the differnt synch info. 


clc 
clear
tic
x=-1; %do not know the factor yet
patientfaktor=[1,-1,-1,1,1,-1,-1,-1,1,1,-1,1,-1,1,-1,x,x,-1]; % for some patients factor 1 or -1 is better. Checked visually
pat=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18];
pat=10;

padding=0; %Determine if the RR should be same length as ECG. Don`t have to be
saving=1;
plotting=0; %plotting Ralphs RR detection

missing_synched_data=[];

Matlabfolder=('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific');
addpath([Matlabfolder '\Synchronizing data']);
addpath('C:\Users\310122653\Documents\PhD\Matlab\R peak detection and HRV\RAlps Rpeak detection')
path='E:\';


for i=1:length(pat)
%% diving into folder structure  
    Neonate=pat(i);
    Ralphsfactor=patientfaktor(1,Neonate);   %Determine if the ECG signal should be turned -1 or not 1. 

    folder=dir([path 'cECG_study\*_test' num2str(pat(i)) ]);     %find the patient folder name...
    datafolder=([path 'cECG_study\' num2str(folder.name) '\']);  % ...go to that folder
    DAQ_patientfolder=dir([datafolder '*_1*']);                  % in the patient folders for DAQ measurements
    
    for k=1:length(DAQ_patientfolder)                            %dive into patient folders 
        ECG_loadfolder=([datafolder DAQ_patientfolder(k,1).name '\']);
        disp('Starting with folder');disp(ECG_loadfolder);
        savepath=[ECG_loadfolder 'Synched Data\' ]; %save the new EDR signals in the same path as the synchronized ECG etc.
        if exist([ECG_loadfolder 'Synched Data\' ], 'dir')==7  % sometimes synch ECG is missing
%% Loading  ECG and creating EDR signal           
            if exist([ECG_loadfolder 'Synched Data\Synced_bin_ECG_500.mat'],'file')
                load([ECG_loadfolder 'Synched Data\Synced_bin_ECG_500.mat']) 
                t_DAQ=linspace(0,length(ECG_bin_synched)/FS_ecg,length(ECG_bin_synched))';
                [DAQ_RR_idx, ~, ~, ~, ~, DAQ_RR, ~] = ecg_find_rpeaks(t_DAQ,Ralphsfactor*ECG_bin_synched, FS_ecg, 250,plotting,0); %, , , maxrate,plotting,saving   -1* because Ralph optimized for a step s slope, we also have steep Q slope. inverting fixes that probel 
                
                % ******* 1# EDR from R peak amplitude not detrended ************

                EDR_FS=length(DAQ_RR_idx)/(length(ECG_bin_synched)/500);                                              % determine sampfle frequency for EDR signal to interpolate to Respiration signal (e.g.: bpm=120 -> sf=2Hz)
                t_edr=linspace(0,floor(length(DAQ_RR_idx)/EDR_FS), length(DAQ_RR_idx))';             % Timeline in seconds with around 2 Hz fs    
                EDR_same_length_as_ECG=interp1(t_edr,DAQ_RR,t_DAQ,'pchip'); 
                EDR_same_length_as_ECG=EDR_same_length_as_ECG-nanmean(EDR_same_length_as_ECG); %center around zero. Remove offset.
                EDR_bin_synched=EDR_same_length_as_ECG; 
                clearvars EDR_same_length_as_ECG
                
                % ***************************************************************
                if saving==1
                    save([savepath 'Synced_bin_EDR_500'],'EDR_bin_synched')
                end
                
            end
            if exist([ECG_loadfolder 'Synched Data\Synced_intelleview_ECG_500Hz.mat'],'file')
                load([ECG_loadfolder 'Synched Data\Synced_intelleview_ECG_500Hz.mat'])
                t_Intellivue=linspace(0,length(ECG_intelleview_synched)/FS_ecg,length(ECG_intelleview_synched))';     
                [Intellivue_RR_idx, ~, ~, ~, ~, Intellivue_RR, ~] = ecg_find_rpeaks(t_Intellivue, Ralphsfactor*ECG_intelleview_synched, FS_ecg, 250,plotting,0); %, , , maxrate,plotting,saving

                % ******* 1# EDR from R peak amplitude not detrended ************

                EDR_FS=length(Intellivue_RR_idx)/(length(ECG_intelleview_synched)/500);                                              % determine sampfle frequency for EDR signal to interpolate to Respiration signal (e.g.: bpm=120 -> sf=2Hz)
                t_edr=linspace(0,floor(length(Intellivue_RR_idx)/EDR_FS), length(Intellivue_RR_idx))';             % Timeline in seconds with around 2 Hz fs    
                EDR_same_length_as_ECG=interp1(t_edr,Intellivue_RR,t_Intellivue,'pchip');
                EDR_same_length_as_ECG=EDR_same_length_as_ECG-nanmean(EDR_same_length_as_ECG);    %center around zero. Remove offset.            
                EDR_intelleview_synched=EDR_same_length_as_ECG;
                clearvars EDR_same_length_as_ECG
                
                % ***************************************************************
                
                if saving==1
                    save([savepath 'Synced_intelleview_EDR_500Hz'],'EDR_intelleview_synched')
                end

            end           

            
         
%% Collect missing folders
        else
            if exist([ECG_loadfolder 'Synched Data\'], 'dir')==0
                missing_synched_data{i,k}=ECG_loadfolder;
            end
            if exist([ECG_loadfolder 'cECGand movement files\' ], 'dir')==0
                missing_cECG_data{i,k}=ECG_loadfolder;
            end
        end % if synch/cECG exist
    disp(['finished with folder ' num2str(ECG_loadfolder)])
    end % for every folder of patient
disp(['finished with patient ' num2str(pat(1,i))])
end % for each patient
if isempty(missing_synched_data)==0
 missing_synched_data=missing_synched_data(~cellfun('isempty',missing_synched_data))  ;
 missing_cECG_data=missing_cECG_data(~cellfun('isempty',missing_cECG_data))  ;
 
 
disp('Finished')
disp(' synch folders are missing for: ')
% for j=1:length(missing_synched_data)
    disp(missing_synched_data)
disp(' cECG folders are missing for: ')
    disp(missing_cECG_data)
else
    disp('finished creating RR signals. ')
    disp('no synched data folder is missing :-)')
end
    toc
    pause(5)
beep