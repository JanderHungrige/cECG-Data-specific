% create cECG and movement files
clear
tic
matlabpath=('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Create data\creating cECG and movement files\');
addpath ('C:\Users\310122653\Documents\PhD\cECG Data\Matlab\creating cECG and movement files\cECG_toolbox')
addpath('C:\Users\310122653\Documents\PhD\cECG Data\Matlab\creating cECG and movement files\cECG_toolbox\cECG_display_functions')
addpath('C:\Users\310122653\Documents\PhD\cECG Data\Matlab\creating cECG and movement files\cECG_toolbox\cECG_processing_functions')
dpath='e:\cECG_study\A_RawData\';


patient_folder_DAQ='20120704_test4';
DAQ_folder='20120704_125601';


cd([dpath patient_folder_DAQ '\'])
      


dataSetpath{1} = ([dpath patient_folder_DAQ '\' DAQ_folder]);
dataSetpath{2} = ([dpath patient_folder_DAQ '\' DAQ_folder '\cECG and movement files\']);
cd(dataSetpath{1})

mkdir('cECG and movement files') % create synched data folder
cd(dataSetpath{2})
mkdir([dataSetpath{2} 'cECG_analysis_results']);                
dataSetpath{3}=([dataSetpath{2} 'cECG_analysis_results']);
cd(dataSetpath{1})


bin_files=dir('*ECG1*.bin'); % find all REf ECG files
numFiles=size(bin_files,1);%amount of ECG ref files

% for <loop over all folders>
dataPath = dataSetpath{1}; % <fill in the folder to be analyzed > (for example: 'D:\cECG\raw_DATA\20120704_test4\20120704_182322')
% numFiles = ; % <fill in the number of bin files in that folder> (for example: 1)
outputPath = (dataSetpath{2}); % <fill in the folder where the resulting signals should be saved> (for example: 'D:\cECG\processed_DATA\20120704_test4\20120704_182322')
cd(matlabpath)
launch_processing

            

% clearvars -except pat dpath patient_folder_DAQ sessions o p    

