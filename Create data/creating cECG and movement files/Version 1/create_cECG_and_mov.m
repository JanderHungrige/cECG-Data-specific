function create_cECG_and_mov()
%SCRIPT for Jan Werth
%
% This script extracts the capacitve ECG + R peak positions + a measure of
% the motion level in the mat files called:
% >> 'motion_level.mat' (250 Hz)
% >> 'cap_ECG.mat' (250 Hz) 
% >> 'cap_Rpeakposition.mat' (in number of samples)
% >> 'ref_Rpeakposition.mat' (in number of samples)
%
% NOTE: ALL the signals above are delayed by 1.74 seconds compared to the
% reference ECG !! (also the ref_Rpeakpositions!!!)
%
%  >> 'referenceECG'(250 Hz)
%
%
% aline.serteyn@gmail.com
% 13-4-2017
tic
clear
addpath ('C:\Users\310122653\Documents\PhD\cECG Data\Matlab\creating cECG and movement files\cECG_toolbox')
addpath('C:\Users\310122653\Documents\PhD\cECG Data\Matlab\creating cECG and movement files\cECG_toolbox\cECG_display_functions')
addpath('C:\Users\310122653\Documents\PhD\cECG Data\Matlab\creating cECG and movement files\cECG_toolbox\cECG_processing_functions')

dpath='e:\cECG_study\';
% patient_foler_DAQ='20120510_test2';
% DAQ_folder='20120510_123507';

cd(dpath)
pat=dir('*_test*');
for o=5:(size(pat,1))-1 %2 because 1 is finsihed, -1 because 15 is different
    patient_foler_DAQ=pat(o).name;
    cd([dpath patient_foler_DAQ '\'])

    sessions=dir ('2012*_1*');
    for p=9:size(sessions,1) % going throug each session of each patient
        if o==4 || o==6 || o==15 % skipping pat 4 6 and 15 as they have errors
            skip=1;
        else 
            skip=0;
        end
        
        if skip==0
             DAQ_folder=sessions(p).name;

            dataSetpath{1} = ([dpath patient_foler_DAQ '\' DAQ_folder]);
            dataSetpath{2} = ([dpath patient_foler_DAQ '\' DAQ_folder '\cECG and movement files\']);
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

            launch_processing
            clearvars -except pat dpath patient_foler_DAQ sessions o p
        end %skip
    %% Display results for one folder
    FS = 500; 
    % figure; aax(1) = subplot(4,1,1); 
    % plot([1:length(cap_ECG)]./FS, cap_ECG); hold on; 
    % plot(cap_Rpeakposition./FS, cap_ECG(cap_Rpeakposition), 'ok'); title('Capacitive ECG and detected peaks')
    % aax(2) = subplot(4,1,2);
    % plot(RRpositions(1:end-1)./FS, diff(RRpositions)./FS,'k'); hold on
    % plot(RRpositionsREF(1:end-1)./FS, diff(RRpositionsREF)./FS,'r'); title('RR tacogram'); ylabel('RR interval [s]'); legend('cap ECG', 'ref ECG'); ylim([0.2 0.9]) %HR 66 tot 300 bpm
    % aax(3) = subplot(4,1,3);
    % plot([1:length(motion_level)]./FS, motion_level); title('Motion level (low motion = below 1e8, high motion = above 1e9)')
    % xlabel('Time [s]')
    % aax(4) = subplot(4,1,4); 
    % ref_RpeakpositionSYNC = max(1,ref_Rpeakposition-FS*1.74);
    % %%% If you have the reference ECG:
    % plot([1:length(referenceECG)]./FS, referenceECG); hold on
    % plot(ref_RpeakpositionSYNC./FS, referenceECG(ref_RpeakpositionSYNC), 'or'); title('Reference ECG and detected peaks')
    % %%% Otherwise:
    % % plot(ref_RpeakpositionSYNC./FS, ones(1,length(ref_RpeakpositionSYNC)), 'or'); title('Reference ECG and detected peaks')
    % xlabel('Time [s]')
    % linkaxes(aax, 'x')
     end%sessions
end % pat
end
toc

