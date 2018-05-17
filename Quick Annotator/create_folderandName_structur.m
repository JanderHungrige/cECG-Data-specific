

function create_folderandName_structur(participant,video)


videoformat='avi'; %change if videoforamt is different

drive='e';
Datapath=[drive ':\cECG_study\A_RawData\'];
savepath=[drive ':\cECG_study\For_Quick_Annotator\new\'];
% participant=12;
cd(Datapath)
folder=dir(['*_test' num2str(participant)]);
cd([Datapath folder.name])

QAfolder='C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Quick Annotator';

%% Create folders for QA output annotations
if exist([savepath 'Annotations'], 'file')==0
status = mkdir ([savepath 'Annotations']);
end
if exist([savepath 'Annotations\Participant' num2str(participant)], 'file')==0
status = mkdir ([savepath 'Annotations\Participant' num2str(participant)]);
end
%% Create folders
%*******************************************************************
%loading the names of the foldes to convert them later to unix code
%*******************************************************************
 names=dir('2012*_1*'); % find all session folders
    numFiles=size(names,1);%amount of folders minus video folder
    if numFiles ~=0 %if there is at least one file
        for fileNumber=1:numFiles % go through each session of particular participant
            if exist([Datapath folder.name '\' names(fileNumber,1).name '\Synched Data'],'file')
                 sessiondate_cell{fileNumber}= names(fileNumber,1).name;  % Read single name in               
                 disp(['Working on Session ' num2str(names(fileNumber,1).name)])
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

                %Create new folder if not existing
                %Go to folder
                if exist([savepath 'participant' num2str(participant)],'file')==0
                    status = mkdir([savepath 'participant' num2str(participant)]);  
                end
                if exist([savepath 'participant' num2str(participant) '\camera\' date],'file')==0
                    status = mkdir([savepath 'participant' num2str(participant) '\camera\' date]);
                end
                if exist([savepath 'participant' num2str(participant) '\DAQ\' date],'file')==0
                    status = mkdir([savepath 'participant' num2str(participant) '\DAQ\' date]);
                end
                if exist([savepath 'participant' num2str(participant) '\Intellivue\' date],'file')==0
                status = mkdir([savepath 'participant' num2str(participant) '\Intellivue\' date]);
                end

                %% Load data - create cell - didstribute to folder          
                sessiondate_cell{fileNumber}= names(fileNumber,1).name;  % Read single folder name in ...

    % load the BIN synched files
    cd([Datapath folder.name])

                cd([Datapath folder.name '\' sessiondate_cell{fileNumber} '\Synched data']) % go to that folder
                if exist([pwd '\Synced_bin_ECG_500.mat'],'file')==2
                    load('Synced_bin_ECG_500'); 
                    load('Synced_bin_EDR_500b');
                    load('Synced_bin_Resp_500');
                    load('Synced_motion_level');
                    load('Synced_cECG_500');
                    load('Synced_bin_RR_500');

        % Create new BIN variables
                    cd(QAfolder)
                    ECG = qa_create_signal('ECG_DAQ', 'uV', 1/FS_ecg, ECG_bin_synched);
                    EDR = qa_create_signal('EDR_DAQ', 'uV', 1/FS_ecg, EDR_bin_synched);
                    Resp =qa_create_signal('Resp_DAQ', 'uV', 1/FS_ecg, Resp_bin_synched);
                    cECG= qa_create_signal('cECG_DAQ', 'uV', 1/FS_ecg, cap_ECG_synched);
                    Motion= qa_create_signal('Motion_DAQ', 'uV', 1/FS_ecg, motion_level_synched);
                    HRV=qa_create_signal('HRV_DAQ', 'uV', 1/FS_ecg, Synced_bin_RR_500);
                    % create filename
                    filename = sprintf('DAQ_%d_%d.mat', Unix1, participant);
                    % Store in new MAT-file
                    save(([savepath 'participant' num2str(participant) '\DAQ\' date '\' filename]), 'ECG', 'EDR','Resp', 'cECG', 'Motion','HRV');
                else
                    disp(['No Bin ECG file found for participant' num2str(participant) ' session ' sessiondate_cell{fileNumber}])                
                end

    % load the IntelliView synched files
                cd([Datapath folder.name '\' sessiondate_cell{fileNumber} '\Synched data']) % go to that folder
                if exist([pwd '\Synced_intelleview_ECG_500Hz.mat'],'file')==2 % pwd + indentify current folder
                    load('Synced_intelleview_ECG_500Hz'); 
                    load('Synced_intelleview_EDR_500bHz');
                    load('Synced_intelleview_Resp_500Hz');
                    load('Synced_Intellivue_RR_500');
    % Create new Intelleview variables           
                    cd('C:\Users\310122653\Documents\PhD\cECG Data\Quick Annotator')
                    ECG = qa_create_signal('ECG_intellivue', 'uV', 1/FS_ecg, ECG_intelleview_synched);
                    EDR = qa_create_signal('EDR_intellivue', 'uV', 1/FS_ecg, EDR_intelleview_synched);
                    Resp =qa_create_signal('Resp_intellivue', 'uV', 1/FS_ecg, Resp_intelleview_synched);
                    HRV=qa_create_signal('HRV_Intellivue', 'uV', 1/FS_ecg, Synced_Intellivue_RR_500);
                    
                    % create filename
                    filename = sprintf('Intellivue_%d_%d.mat', Unix1, participant);
                    % Store in new MAT-file
                    save(([savepath 'participant' num2str(participant) '\Intellivue\' date '\' filename]), 'ECG', 'EDR','Resp','HRV');
                 else
                    disp(['No intellivue ECG file found for participant' num2str(participant) ' session ' sessiondate_cell{fileNumber}])
                end % if exist intellivue
            end % if exist Synched Data
            cd([Datapath folder.name])
        end
    end
%% Renaming and shifting Video Files
if video==1
cd([Datapath folder.name])
names=dir('*video*');
cd([Datapath folder.name '\' names.name]) %go into videofolder
names=dir(['*.' videoformat]);

    numFiles=size(names,1);%amount of folders minus video folder
    if numFiles ~=0 %if there is at least one file
        for fileNumber=1:numFiles % go through each session of particular participant
            sessiondate_cell{fileNumber}= names(fileNumber,1).name;  % Read single name in         
            str=sessiondate_cell{1,fileNumber};
            ss=strsplit(str,'.');
            s = strsplit(str,'_'); 

            dateStr = s{1}; 
            year = dateStr(1:4); 
            month = dateStr(5:6); 
            day = dateStr(7:8); 
            date = [year '-' month '-' day]; 
            % Conert time into Unix time
            Datum = datetime(ss{1}, 'InputFormat', 'yyyyMMdd_HHmmss');    
            Datum.TimeZone = 'Europe/Berlin';
            Unix1 = posixtime(Datum);
            %rename 
            copyfile( names(fileNumber).name,[str '#backup'] ) %create backup of video with old name
            movefile( names(fileNumber).name, sprintf(['%i.%i.' videoformat],participant, Unix1));%change the filenames
            newnames=dir([num2str(participant) '.*']);
            %Moving
            movefile( newnames.name ,[savepath 'participant' num2str(participant) '\camera\' date])
        end
    end
    
 %% rename backupvideos back to old name  
cd([Datapath folder.name ])
names=dir('*video*');
cd([Datapath folder.name '\' names.name]) %go into videofolder
names=dir('*backup*');
    numFiles=size(names,1);%amount of folders minus video folder
    if numFiles ~=0 %if there is at least one file
        for fileNumber=1:numFiles % go through each session of particular participant
            sessiondate_cell{fileNumber}= names(fileNumber,1).name;  % Read single name in         
            str=sessiondate_cell{1,fileNumber};
            s = strsplit(str,'#'); 

            oldname = s{1}; 
            %rename 
            movefile(names(fileNumber).name, sprintf('%s',oldname));%change the filenames
        end
    end
end %if video

end




