%Create RR signal from DAQ and/or Intellivue signal

clc 
clear
tic

pat=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18];
Ralphsfactor={1;-1;-1;1;1;-1;-1;-1;1;1;-1;1;-1;1;-1;-1;-1;-1};%Determine if the ECG signal should be turned -1 or not 1. 
%***************************************************************************
% Attention :
% PATIENT 16 HAS ONEFOLDER STARTING WITH _09 INSTEAD OF _1. cHANGE THAT MANUALLY FOR ONE RUN
%***************************************************************************

padding=0; %Determine if the RR should be same length as ECG. Don`t have to be
saving=1;
plotting=1; %plotting Ralphs RR detection

missing_synched_data=[];

Matlabfolder=('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific');
addpath([Matlabfolder '\Synchronizing data']);
addpath('C:\Users\310122653\Documents\PhD\Matlab\R peak detection and HRV\RAlps Rpeak detection')
path='E:\';


for i=1:length(pat)
%% diving into folder structure  
    Neonate=pat(i);
    folder=dir([path 'cECG_study\*_test' num2str(pat(i)) ]);     %find the patient folder name...
    datafolder=([path 'cECG_study\' num2str(folder.name) '\']);  % ...go to that folder
    DAQ_patientfolder=dir([datafolder '*_1*']);                  % in the patient folders for DAQ measurements
    
    for k=1:length(DAQ_patientfolder)                            %dive into patient folders 
        ECG_loadfolder=([datafolder DAQ_patientfolder(k,1).name '\']);
        disp('Starting with folder');disp(ECG_loadfolder);
        savepath=[ECG_loadfolder 'Synched Data\' ]; %save the new RR signals in the same path as the synchronized ECG etc.
        if exist([ECG_loadfolder 'Synched Data\' ], 'dir')==7  % sometimes synch ECG is missing
%% Loading  ECG and creating RR signal           
            if exist([ECG_loadfolder 'Synched Data\Synced_bin_ECG_500.mat'],'file')
                load([ECG_loadfolder 'Synched Data\Synced_bin_ECG_500.mat']) 
                t_DAQ=linspace(0,length(ECG_bin_synched)/FS_ecg,length(ECG_bin_synched))';
                Synced_bin_RR_500=nan(1,length(ECG_bin_synched));
                [DAQ_RR_idx, ~, ~, ~, ~, DAQ_RR, ~] = ecg_find_rpeaks(t_DAQ,Ralphsfactor{Neonate,1}*ECG_bin_synched, FS_ecg, 250,plotting,0); %, , , maxrate,plotting,saving   -1* because Ralph optimized for a step s slope, we also have steep Q slope. inverting fixes that probel 
                figure;title('DAQ')
                Synced_bin_RR_500(DAQ_RR_idx)=DAQ_RR;
                if saving==1
                    save([savepath 'Synced_bin_RR_500'],'Synced_bin_RR_500')
                end
            end
            if exist([ECG_loadfolder 'Synched Data\Synced_intelleview_ECG_500Hz.mat'],'file')
                load([ECG_loadfolder 'Synched Data\Synced_intelleview_ECG_500Hz.mat'])
                t_Intellivue=linspace(0,length(ECG_intelleview_synched)/FS_ecg,length(ECG_intelleview_synched))';     
                Synced_Intellivue_RR_500=nan(1,length(ECG_intelleview_synched));                
                [Intellivue_RR_idx, ~, ~, ~, ~, Intellivue_RR, ~] = ecg_find_rpeaks(t_Intellivue, Ralphsfactor{Neonate,1}*ECG_intelleview_synched, FS_ecg, 250,plotting,0); %, , , maxrate,plotting,saving
                figure;title('Intellivue')                
                Synced_Intellivue_RR_500(Intellivue_RR_idx)=Intellivue_RR;
                if saving==1
                    save([savepath 'Synced_Intellivue_RR_500'],'Synced_Intellivue_RR_500')
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
beep-