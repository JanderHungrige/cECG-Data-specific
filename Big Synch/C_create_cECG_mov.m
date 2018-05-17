% create cECG and movement files

% This m FIle uses Version 1 to load the cECG and movement files. This does
% not work for patients with more than 1 ECG file. Better try Version 2

% in C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Create data\creating cECG and movement files\Version 1
% for the manual oneby one file I repaired the read_data to be able to read
% patients with more than one ECG file

tic
clear
matlabpath=('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Create data\creating cECG and movement files\');

% addpath ('C:\Users\310122653\Documents\PhD\cECG Data\Matlab\creating cECG and movement files\cECG_toolbox')
% addpath('C:\Users\310122653\Documents\PhD\cECG Data\Matlab\creating cECG and movement files\cECG_toolbox\cECG_display_functions')
% addpath('C:\Users\310122653\Documents\PhD\cECG Data\Matlab\creating cECG and movement files\cECG_toolbox\cECG_processing_functions')

addpath('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Create data\creating cECG and movement files\Version 1')
addpath('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Create data\creating cECG and movement files\Version 1\cECG_processing_functions')
addpath('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Create data\creating cECG and movement files\Version 1\cECG_display_functions')

dpath='e:\cECG_study\';


SAVE=0; % determines if the file offline_analysis from ALine is started. There the cECG is created from the 8 sensors. IF you only need the Raw cECG that is not neccesary.
saving_RAW_cECG=1; % Used in read_data from Serteyn, Aline <aline.serteyn@philips.com>; Use that to save the RAW cECG for the algorythm from Guy Warmerdam (G.J.J.Warmerdam@tue.nl) 

% patient_foler_DAQ='20120905_test15';
% DAQ_folder='20120905_121501';

cd(dpath)
pat=dir('*_test*');
for o=4%(size(pat,1)) %2 because 1 is finsihed, -1 because 15 is different
    disp(['working on patient' num2str(pat(o,1).name)])
    patient_folder_DAQ=pat(o).name;
    cd([dpath patient_folder_DAQ '\'])

    sessions=dir ('2012*_1*');
    for p=1:size(sessions,1) % going throug each session of each patient
        disp(['working on session' num2str(sessions(p,1).name)])

        if o==4 || o==6 || o==15 || o==16|| o==17|| o==18  % skipping pat 4 6 and 15 as they have errors
            skip=0;
        else 
            skip=0;
        end
        
        if o~=9
            skip=0;
        elseif o==9
            skip=1;
        end
        if skip==0
             DAQ_folder=sessions(p).name;

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

            
        end
        clearvars -except pat dpath patient_folder_DAQ sessions o p    
    end%sessions
end % pat