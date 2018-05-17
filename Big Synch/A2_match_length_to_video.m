% USE THIS AFTER YOU CREATED AND AYNCHED ALSO cECG AND MOV WITH 
%C_synching_cECG_to_data

%Make the ECG and Resp at least the same length as the video. Data can be
% longer. This is important for the later spplitting into 1 hour. If the
% data is to short, an missmatching error is thrown. 
clear
clc
 pat=[1,2,3,4];
 
 
Matlabfolder=('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific');
path='E:\';
folderpath='cECG_study\';
folder=dir([path 'cECG_study\*_test*' ]);


for i=1:length(pat)
    patientfolder=dir([path 'cECG_study\*_test' num2str(pat(1,i)) ]);
    cd([path folderpath patientfolder.name])
    sessionfolder=dir('2012*_1*');
    for k=1:length(sessionfolder)
       cd([path folderpath patientfolder.name '\' sessionfolder(1,k).name '\Synched Data' ])  
       variablenames=dir('*.mat');
    
        for j=1:length(variablenames)
            load(variablenames(j,1).name)
        end
     
        % while beeing in the session folder also load the video
        % information
    cd([path folderpath patientfolder.name])
    videofolder=dir('*_video_*');
    cd([path folderpath patientfolder.name '\' videofolder.name]) 
    vidInfo=VideoReader([sessionfolder(1,k).name '.avi']); 
    
    %compare the duration of video and ECG
    if length(ECG_bin_synched)/500 < vidInfo.Duration % if Video is longer than ECG...
        addSamples=nan(ceil((vidInfo.Duration-length(ECG_bin_synched)/500)*500),1); %ad this amount of samples  
        for j=1:length(variablenames)
            changevars= who( '*_synched');
            eval(['vertcat(',changevars{j,1}, addSamples ')']); 
        end
    
    end
    
    end  % loading sessions

    
    
end % for each patient


