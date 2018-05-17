

clc 
clear
clc

pat=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18];
% pat=[1,2,3,5,7,8,9,10,11,12,13,14,15];
pat=4;

padding=0; %Determine if the cECG and mov should be same length as ECG. Don`t have to be

missing_synched_data=[];

Matlabfolder=('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific');
synchfolder=([Matlabfolder '\Synchronizing data']);
path='E:\';

folder=dir([path 'cECG_study\A_RawData\*_test*']);



for i=1:length(pat)
  
    Neonate=pat(i);
    datafolder=([path 'cECG_study\A_RawData\' num2str(folder(Neonate,1).name) '\']); % All patient folders
    DAQ_patientfolder=dir([datafolder '*_1*']);                      % in the patient folders for DAQ measurements
    intelliview_patientfolder=dir([datafolder '*Intelliview*']);     % in the patient folders for Intelliview measurements
    
    for k=1:length(DAQ_patientfolder)
        
        DAQ_loadfolder=([datafolder DAQ_patientfolder(k,1).name '\']);
%         if isempty(intelliview_patientfolder)==0 % Patient 1,2 and 3 do not have intelliview data
%         intelliview_loadfolder=intelliview_patientfolder(k,1).name;
%         else
%         intelliview_loadfolder=[];
%         end

        if exist([DAQ_loadfolder 'Synched Data\' ], 'dir')==7 & exist([DAQ_loadfolder 'cECG and movement files\' ], 'dir')==7 % sometimes synch ECG is missing, sometime the cECG
            cd(synchfolder)
                Synch_cECG_and_movement_to_ECG(DAQ_loadfolder)
                
%make the cECG and mov the same length as the ECG
            if padding==1
                load([DAQ_loadfolder 'Synched Data\Synced_bin_ECG_500']);
                load([DAQ_loadfolder 'Synched Data\Synced_cECG_500']);
                load([DAQ_loadfolder 'Synched Data\Synced_motion_level']);
                
                if length(ECG_bin_synched)>length(cap_ECG_synched)
                    length(ECG_bin_synched)-length(cap_ECG_synched);
                    matching=nan(length(ECG_bin_synched)-length(cap_ECG_synched),1);
                    cap_ECG_synched=[cap_ECG_synched; matching]; % adding nans at the end to gain same length
                    save([DAQ_loadfolder 'Synched Data\Synced_cECG_500'],'cap_ECG_synched');
                elseif length(ECG_bin_synched)<length(cap_ECG_synched)
                    cap_ECG_synched=cap_ECG_synched(1:length(ECG_bin_synched),1);
                    save([DAQ_loadfolder 'Synched Data\Synced_cECG_500'],'cap_ECG_synched');                    
                end
                    
                if length(ECG_bin_synched)>length(motion_level_synched)
                    length(ECG_bin_synched)-length(motion_level_synched);
                    matching=nan(length(ECG_bin_synched)-length(motion_level_synched),1);
                    motion_level_synched=[motion_level_synched; matching]; % adding nans at the end to gain same length
                    save([DAQ_loadfolder 'Synched Data\Synced_motion_level'],'motion_level_synched');                  
                elseif length(ECG_bin_synched)<length(motion_level_synched)
                    motion_level_synched=motion_level_synched(1:lenght(ECG_bin_synched),1);
                    save([DAQ_loadfolder 'Synched Data\Synced_motion_level'],'motion_level_synched');              
                end  
            end
% Collect missing folders
        else
            if exist([DAQ_loadfolder 'Synched Data\'], 'dir')==0
                missing_synched_data{i,k}=DAQ_loadfolder;
            end
            if exist([DAQ_loadfolder 'cECGand movement files\' ], 'dir')==0
                missing_cECG_data{i,k}=DAQ_loadfolder;
            end
        end % if synch/cECG exist
    disp(['finished with folder ' num2str(DAQ_loadfolder)])
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
    disp(['finished synching and extending/cutting cECG and movement to ECG. '])
    disp('no synched data folder is missing :-)')
end
    