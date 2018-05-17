clear 
clc
tic

Matlabfolder=('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific');
addpath([Matlabfolder '\Create data']);
addpath('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\R peak detection and HRV\RAlps Rpeak detection');
path='E:\';
folder=dir([path 'cECG_study\*_test*' ]);

% ***********************************************************************
pat=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18];
pat=[8,9,10,11,12,13,14,15,16,17,18];
pat=[9]

x=-1; %do not know the factor yet
patientfaktor=[1,-1,-1,1,1,-1,-1,-1,1,1,-1,1,-1,1,-1,x,x,-1]; % for some patients factor 1 or -1 is better. Checked visually
saveing=1;
% ***********************************************************************

%example
% patient_foler='20120724_test10';
% DAQ_folder='20120724_112103';
% Intellivue_folder='2012-07-24 11-21_001';



for i=1:length(pat)
  
    Neonate=pat(i);
    disp('***************************************************************')
    disp(['Working on Pat:' num2str(Neonate)])
    factor=patientfaktor(Neonate);
    disp(['using factor: ' num2str(factor)]);
    patient_folder=([path 'cECG_study\' num2str(folder(Neonate,1).name) '\']); % All patient folders
    DAQ_folders=dir([patient_folder '*_1*']);                      % in the patient folders for DAQ measurements
    intellivue_patientfolder=dir([patient_folder '*Intelliview*']);     % in the patient folders for Intelliview measurements
    Intellivue_folders=dir([patient_folder '\' intellivue_patientfolder.name '\2012*']);
    
    for k=1:length(DAQ_folders)
        DAQ_folder=DAQ_folders(k,1).name;
        disp(['* working on ' num2str(DAQ_folder)])
        dataSetpath{1}=([patient_folder DAQ_folder '\']);
        
        if isempty(Intellivue_folders)==0 % Patient 1,2 and 3 do not have intelliview data
            intellivue=1;
            intellivue_folder=Intellivue_folders(k,1).name;
            
            dataSetpath{2}=([patient_folder intellivue_patientfolder.name '\' intellivue_folder '\']);
        else
            intellivue=0;
        end
        Respiration_from_ECG_from_original(1,intellivue,saveing,factor,dataSetpath);
        
    end
               
                
                
end
toc