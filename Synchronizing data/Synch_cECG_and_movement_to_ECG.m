% Synch cECG and movement to the reference ECG 
% 
% Note that the capacitive ECG is delayed by 1.74 seconds in comparison with the reference ECG.
% This file loads the *.mat file of the synched reff file and the info on
% the synch. 
% Then data is added accordantly to the synchronisation of the ref ECG and
% the 1.74s delay. 
function Synch_cECG_and_movement_to_ECG(loadfolder)

% folder=dir([loadfolder '*_test']);
% dpath='I:\cECG_study\';
% patient_foler_DAQ='20120504_test1';
% DAQ_folder='20120504_124005';



% dataSetpath{1} = ([dpath patient_foler_DAQ '\' DAQ_folder '\Synched Data']);
% cd(dataSetpath{1})
% dataSetpath{2} = ([dpath patient_foler_DAQ '\' DAQ_folder '\cECG and movement files\']);
% addpath(dataSetpath{2})
dataSetpath{1}=([loadfolder 'Synched Data']);
dataSetpath{2}=([loadfolder 'cECG and movement files\']);
cd(dataSetpath{1})
addpath(dataSetpath{2})

load('Synchronization_info', '-mat','offsetCAM_500')
load([dataSetpath{2} 'cap_ECG'], '-mat')
load([dataSetpath{2} 'cap_Rpeakposition'], '-mat')
load([dataSetpath{2} 'motion_level'], '-mat')
load([dataSetpath{2} 'ref_Rpeakposition'], '-mat')



delay=1.74*500;
newoffset=offsetCAM_500+delay;

if newoffset >0 %Camera starts later then ECG
    %add nan or zeroes to cECG
    addition(1:newoffset,1)=NaN;
    cap_ECG_synched        =[addition' cap_ECG];
    motion_level_synched   =[addition'  motion_level] ;
    cap_Rpeakposition      =cap_Rpeakposition+newoffset;
elseif newoffset <0 %camera started earlier than ECG
    %cut the first few samples of ECG, intelleview ECG, Respiration
    cap_ECG_synched         =cap_ECG(-newoffset:end,1);
    motion_level_synched    =motion_level(-newoffset:end,1);
    cap_Rpeakposition       =cap_Rpeakposition+newoffset;
elseif newoffset==0
    addition(1:delay,1)=NaN;
    cap_ECG_synched        =cap_ECG;
    motion_level_synched   =motion_level;
    cap_Rpeakposition      =cap_Rpeakposition;
    disp('Cam and ECG are same Please check if data is Ok')
end

cap_ECG_synched=cap_ECG_synched'; motion_level_synched=motion_level_synched';

 savedatapath=([dataSetpath{1}]);
 save([savedatapath '\Synced_cECG_500'],'cap_ECG_synched')
 save([savedatapath '\Synced_motion_level'],'motion_level_synched')
 save([savedatapath '\Synced_cap_Rpeakposition'],'cap_Rpeakposition')
end
