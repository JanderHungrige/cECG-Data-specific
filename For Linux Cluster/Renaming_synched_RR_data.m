% Renaming single files to new time format
% Unfortunattely, some file where not included in the files for the quicl
% annotator. Eg.the RR peaks from the DAQ system. Therefore this RR files
% have to be renamed separaetluy to know to which session they belong.

clear
clc
participant=4;


for P=1:length(participant)
    drive='E:\';
    cd([drive 'cECG_study\'])
    folder=dir(['*_test' num2str(participant(P))]);
    cd([drive 'cECG_study\' folder.name])

    savefolder=['E:\cECG_study\Annotations\Datafiles\For Quick Annotator\participant' num2str(participant(P)) '\'];

    % QAfolder='C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Quick Annotator';

    %% Create folders for QA output annotations
    if exist([drive 'cECG_study\For_Quick_Annotator\Annotations'], 'file')==0
    status = mkdir ([drive 'cECG_study\For_Quick_Annotator\Annotations']);
    end
    if exist([drive 'cECG_study\For_Quick_Annotator\Annotations\Participant' num2str(participant(P))], 'file')==0
    status = mkdir ([drive 'cECG_study\For_Quick_Annotator\Annotations\Participant' num2str(participant(P))]);
    end
    %% Create folders
    %*******************************************************************
    %loading the names of the foldes to convert them later to unix code
    %*******************************************************************
     names=dir('2012*_1*'); % find all session folders
        numFiles=size(names,1);%amount of folders minus video folder
        if numFiles ~=0 %if there is at least one file
            for fileNumber=1:numFiles % CHANGED  1:numFiles %go through each session of particular participant
                if exist([drive 'cECG_study\' folder.name '\' names(fileNumber,1).name '\Synched Data'],'file')
                     sessiondate_cell{fileNumber}= names(fileNumber,1).name;  % Read single name in               
                    % Convert name into date for folders
                    str=sessiondate_cell{1,fileNumber};
                    s = strsplit(str,'_'); 

                    dateStr = s{1}; 
                    year = dateStr(1:4); 
                    month = dateStr(5:6); 
                    day = dateStr(7:8); 
                    date = [year '-' month '-' day]; 
                    % Conert time into Unix time
                    Datum = datetime(str, 'InputFormat', 'yyyyMMdd_HHmmss');    
                    Datum.TimeZone = 'Europe/Berlin';
                    Unix1 = posixtime(Datum);


                    % Load / save file
                    load([drive 'cECG_study\' folder.name '\' names(fileNumber,1).name '\Synched Data\Synced_bin_RR_500'])
                    filename = sprintf('DAQ_%d_%d_RR.mat', Unix1,participant(P));
                    datefolder=dir([savefolder 'DAQ\2012*']);
                    save([savefolder 'DAQ\' datefolder.name '\' filename] , 'Synced_bin_RR_500');
                    
                    if exist([drive 'cECG_study\' folder.name '\' names(fileNumber,1).name '\Synched Data\Synced_Intellivue_RR_500.mat'])
                        load([drive 'cECG_study\' folder.name '\' names(fileNumber,1).name '\Synched Data\Synced_Intellivue_RR_500'])
                        filename = sprintf('Intellivue_%d_%d_RR.mat', Unix1,participant(P));
                        datefolder=dir([savefolder 'Intellivue\2012*']);
                        save([savefolder 'Intellivue\' datefolder.name '\' filename], 'Synced_Intellivue_RR_500');                    
                    else
                        disp(['no Intellivue in folder ' sessiondate_cell{fileNumber}])
                    end
                    
                    if exist([drive 'cECG_study\' folder.name '\' names(fileNumber,1).name '\Synched Data\Synced_cap_Rpeakposition'])
                        load([drive 'cECG_study\' folder.name '\' names(fileNumber,1).name '\Synched Data\Synced_cap_Rpeakposition'])
                        filename = sprintf('%d_%d_cRR.mat', Unix1,participant(P));
                        datefolder=dir([savefolder 'DAQ\2012*']);
                        save([savefolder 'DAQ\' datefolder.name '\' filename], 'Synced_cap_Rpeakposition');                       
                    else 
                        disp(['no cECG-RR in folder ' sessiondate_cell{fileNumber} ])
                    end

                end
            end
        end
end
                    