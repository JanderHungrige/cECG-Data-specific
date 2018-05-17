
clear 
clc
tic
% in C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Create data\creating cECG and movement files\Version 1
% for the manual oneby one file I repaired the read_data to be able to read
% patients with more than one ECG file

Matlabfolder=('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Create data\creating cECG and movement files\Version 2\');
addpath([Matlabfolder '\Create data']);
addpath('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\R peak detection and HRV\RAlps Rpeak detection');
path='E:\';
folder=dir([path 'cECG_study\*_test*' ]);

% ***********************************************************************
pat=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18];
pat=[8,9,10,11,12,13,14,15,16,17,18];
pat=[9]
% ***********************************************************************


for i=1:length(pat)
  
    Neonate=pat(i);
    disp('***************************************************************')
    disp(['Working on Pat:' num2str(Neonate)])
%     factor=patientfaktor(Neonate);
%     disp(['using factor: ' num2str(factor)]);
   
    patient_folder=([path 'cECG_study\' num2str(folder(Neonate,1).name) '\']); % All patient folders
    % Pat 16 has one folder starting with _09. Repeate once manualy 
    DAQ_folders=dir([patient_folder '*_1*']);  % in the patient folders for DAQ measurements  
     
    for k=1:length(DAQ_folders)
        DAQ_folder=DAQ_folders(k,1).name;
        disp(['* working on ' num2str(DAQ_folder)])
        dataSetpath{1}=([patient_folder DAQ_folder '\']);
        dataSetpath{2} = ([dataSetpath{1} 'cECG and movement files\']);

% ************** Read capacitive sensors ******************
cd(dataSetpath{1})

        bin_files=dir('Sensor 1_*.bin'); % find all REf ECG files
        numFiles=size(bin_files,1);%amount of ECG ref files
        
        for sens=1:8 % initialize
            eval(['SEN' num2str(sens) '= [];'])
        end
        for fileNumber=0:numFiles-1
            for sens=1:8
                m = eval(['memmapfile(''Sensor ',num2str(sens),'_0000000',num2str(fileNumber),'.bin'',''format'',''uint8'');']);
                x = m.data;%(3*(sampleCounter-offsetECG)+1:3*(sampleCounter-offsetECG)+3*numberNewSamples);
                eval(['SEN' num2str(sens) '_tmp = double(typecast(x(1:3:length(x)),''int8'')'')*2^16+double(x(2:3:length(x))'')*256+double(x(3:3:length(x))'');']);
                eval(['SEN' num2str(sens) '= [SEN' num2str(sens) ' SEN' num2str(sens) '_tmp];']) %creating one long cECG file
                eval(['clearvars SEN' num2str(sens) '_tmp']) 
            end
        end
         Raw_cECGMatrix=[SEN1;SEN2;SEN3;SEN4;SEN5;SEN6;SEN7;SEN8]; % Merge to Matrix
% ************** Creating Data ******************
cd(Matlabfolder)        
        
        [cECG1, cECG2, HR, motionA, motionL, Nr_channels_used, raw_cECG, cap_Resp] = cECG_extraction(Raw_cECGMatrix);
            cap_ECG=cECG1;
            cap_ECG_II=cECG2;
            motion_level=motionA;
            
% ************** Saving ******************
         save([dataSetpath{2} 'cap_ECG'],'cap_ECG')
         save([dataSetpath{2} 'cap_ECG_II'],'cap_ECG_II')
         save([dataSetpath{2} 'motion_level'],'motion_level')
         save([dataSetpath{2} 'motionL'],'motionL')
         save([dataSetpath{2} 'Nr_channels_used'],'Nr_channels_used')
         save([dataSetpath{2} 'raw_cECG'],'raw_cECG')
         save([dataSetpath{2} 'cap_Resp'],'cap_Resp')
 
    end
end
cd('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific')% going back