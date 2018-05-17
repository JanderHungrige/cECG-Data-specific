function split_data_in_1hour(participant)
tic

addpath('C:\Users\310122653\Documents\PhD\cECG Data\Quick Annotator')
drive='e';
% participant=1;
path=([drive ':\cECG_study\For_Quick_Annotator\participant' num2str(participant) ]);
dospath=([drive ':/cECG_study/For_Quick_Annotator/participant' num2str(participant) ]); % for ffmpeg

camPath=([path '\camera\']);
doscamPath=([dospath '/camera/']); %for ffmpeg
DAQPath=([path '\DAQ\']);
IntellivuePath=([path '\Intellivue\']);

addpath(camPath)
cd(camPath)

folder=dir('2012*');


numFiles=size(folder,1);
for i=1:numFiles
%GATHER INFORMATION
    cd([camPath folder(i,1).name])
    status = mkdir([camPath folder(i,1).name], 'splittet_data');                
    names=dir('*.avi');
    numVideos=size(names,1);
    for j=1:numVideos      
        info = mmfileinfo(names(j,1).name); %Getting length of Video
        durationofVideo=info.Duration/3600; % Get how many parts should be created
%IF LONGER THAN ONE HOUR SPLIT OTHERWISE DO NOTHING
        if info.Duration >= 3600 %longer than 1 hour
            split=3600;
            if durationofVideo-floor(durationofVideo)<0.25 %If the duration of last video would be shorter than 15min use differnt split
                durationofVideo=info.Duration/2700; %split on 45min instead of 1 hour
                split=2700; %indication for name change to take 45min(2700s) instead of 1 hour (3600s)
            end           
%% SPLITTING VIDEO         
            s = strsplit(info.Filename,'.'); 
            for k=0:ceil(durationofVideo)-1 %-1 as for the first run we do not add anything to the beginning
                if split == 3600
                    add_start=k*3600;  
                    add_end=k*3600+3600;
                elseif split ==2700
                    add_start=k*2700;
                    add_end=k*2700+2700;
                end
                
                t_start = int64(str2num(s{1,2})+add_start);
                t_end   = int64(str2num(s{1,2})+add_end);
%                 posix_time  = int64(s{1.2}); 
              
                % Input file
                file_in = sprintf([doscamPath folder(i,1).name '/%u.%u.avi'], ...
                  participant, int64(str2num(s{1,2})));
                file_out = sprintf([doscamPath folder(i,1).name '/splittet_data/%u.%u.avi'], ...
                  participant, int64(str2num(s{1,2})+ add_start));
              
                % Assemble command
                % NOTE: ffmpeg.exe should be located in the work folder
                cmd = sprintf('"C:/ffmpeg/bin/ffmpeg.exe" -i %s -ss %u -t %u -c copy %s',...
                  file_in, add_start, add_end, file_out); 
              
                % Execute command
                status = dos(cmd);
                if (status ~= 0)
                  error('Splitting failed');
                end
            end% splitting the video

%% SPLITTING DATA
%splitting DAQ
            if exist(DAQPath,'file')== 7 % if folder exist
                cd([DAQPath folder(i,1).name])
                status = mkdir([DAQPath folder(i,1).name] , '\splittet_data');                                
                load(['DAQ_' s{1,2} '_' num2str(participant)])
                for m=0:ceil(durationofVideo)-1 % cut data same as Video
                    h=m+1; %Cell indice needs to be at least 1
                if split == 3600
                    add_start=m*3600;  
                    add_end=m*3600+3600;
                elseif split ==2700
                    add_start=m*2700;
                    add_end=m*2700+2700;
                end
                    if m==0
                        ECGsplit{h,1}=ECG.values(1:add_end*500);
                        EDRsplit{h,1}=EDR.values(1:add_end*500 );
                        Respsplit{h,1}=Resp.values(1:add_end*500 );
                        if exist('cECG','var')==1
                        cECGsplit{h,1}=cECG.values(1:add_end*500 );
                        Motionsplit{h,1}=Motion.values(1:add_end*500 );  
                        end
                    elseif m~=0 && h<ceil(durationofVideo)              
                        ECGsplit{h,1}=ECG.values(add_start*500:add_end*500 );
                        EDRsplit{h,1}=EDR.values(add_start*500:add_end*500 );
                        Respsplit{h,1}=Resp.values(add_start*500:add_end*500 );
                        if exist('cECG','var')==1
                        cECGsplit{h,1}=cECG.values(add_start*500:add_end*500 );
                        Motionsplit{h,1}=Motion.values(add_start*500:add_end*500 );
                        end
                    else % if the last file is shorter than 2700 or 3600 
                        ECGsplit{h,1}=ECG.values(add_start*500:end);
                        EDRsplit{h,1}=EDR.values(add_start*500:end);
                        Respsplit{h,1}=Resp.values(add_start*500:end) ;  
                        if exist('cECG', 'var')==1
                        cECGsplit{h,1}=cECG.values(add_start*500:end);
                        Motionsplit{h,1}=Motion.values(add_start*500:end );                        
                        end
                    end
                ECG = qa_create_signal('ECG_DAQ', 'uV', 1/500, ECGsplit{h,1});
                EDR = qa_create_signal('EDR_DAQ', 'uV', 1/500, EDRsplit{h,1});
                Resp =qa_create_signal('Resp_DAQ', 'uV', 1/500, Respsplit{h,1});
                if exist('cECG','var')==1
                cECG= qa_create_signal('cECG_DAQ', 'uV', 1/500, cECGsplit{h,1});
                Motion= qa_create_signal('Motion_DAQ', 'uV', 1/500, Motionsplit{h,1});       
                end
                t_start = int64(str2num(s{1,2})+add_start);
                t_end   = int64(str2num(s{1,2})+add_end);

                % create filename
                filename = sprintf('DAQ_%d_%d.mat', t_start, participant);
                % Store in new MAT-file
                if exist('cECG','var')==1
                save(([ DAQPath folder(i,1).name '\splittet_data\' filename]), 'ECG', 'EDR','Resp','cECG','Motion');
                else
                save(([ DAQPath folder(i,1).name '\splittet_data\' filename]), 'ECG', 'EDR','Resp');
                end
                    
                end % for count of videos to run throug DAQ             
            end % if exist DAQ
            
%% SPLITTING INTELLIVUE

            cd([IntellivuePath folder(i,1).name])
            intelfolder=dir([IntellivuePath 'Intellivue*']);

            if isempty(intelfolder)==0 %check if intellivue files exist
                status = mkdir([IntellivuePath folder(i,1).name], '\splittet_data');                                
                load(['Intellivue_' s{1,2} '_' num2str(participant)])
                for m=0:ceil(durationofVideo)-1 % cut data same as Video
                    h=m+1;%Cell indice needs to be at least 1
                if split == 3600
                    add_start=m*3600;  
                    add_end=m*3600+1*3600;
                elseif split ==2700
                    add_start=m*2700;
                    add_end=m*2700+1*2700;
                end
                    if m==0
                        ECGsplit{h,1}=ECG.values(1:add_end*500);
                        EDRsplit{h,1}=EDR.values(1:add_end*500 );
                        Respsplit{h,1}=Resp.values(1:add_end*500 );                      
                    elseif m~=0 && h<ceil(durationofVideo)                
                        ECGsplit{h,1}=ECG.values(add_start*500:add_end*500 );
                        EDRsplit{h,1}=EDR.values(add_start*500:add_end*500 );
                        Respsplit{h,1}=Resp.values(add_start*500:add_end*500 );

                    else % if the last file is shorter than 2700 or 3600 
                        ECGsplit{h,1}=ECG.values(add_start*500:end);
                        EDRsplit{h,1}=EDR.values(add_start*500:end);
                        Respsplit{h,1}=Resp.values(add_start*500:end) ;      
                      
                    end
                ECG = qa_create_signal('ECG_DAQ', 'uV', 1/500, ECGsplit{h,1});
                EDR = qa_create_signal('EDR_DAQ', 'uV', 1/500, EDRsplit{h,1});
                Resp =qa_create_signal('Resp_DAQ', 'uV', 1/500, Respsplit{h,1});
      
                
                t_start = int64(str2num(s{1,2})+add_start);
                t_end   = int64(str2num(s{1,2})+add_end);

                % create filename
                filename = sprintf('Intellivue_%d_%d.mat', t_start, participant);
                % Store in new MAT-file
                save(([IntellivuePath folder(i,1).name '\splittet_data\' filename]), 'ECG', 'EDR','Resp');

            
            
                end % for count of videos to run throug Intellivue             
            else
                disp('No intellivue files found')
            end % if exist intellivue            
            
            
        else % if duration shorter than 1 hour
            disp(['No splitting neccecary. File ' names(j,1).name ' shorter than 1 hour'])
        end % if longer than 1 hour
            cd([camPath folder(i,1).name]) %back to main 

    end % for several videos
%% SWITHCING SPLITTET AND LONG FILES
% CREATING BACKUP FOLDERS 
%     status = mkdir([camPath folder(i,1).name] , '\Long_data');                                
%     status = mkdir([DAQPath folder(i,1).name] , '\Long_data');                                
%     status = mkdir([IntellivuePath folder(i,1).name] , '\Long_data');
    
% Copy original files into backupfolder
    if isempty(dir([camPath folder(i,1).name '\*.avi']))==0    
        copyfile ([camPath folder(i,1).name '\*.avi'] ,[camPath folder(i,1).name '\Long_data'])
    else
        disp('No Videos found')
    end
    if isempty(dir([DAQPath folder(i,1).name '\*.mat']))==0        
        copyfile ([DAQPath folder(i,1).name '\*.mat'],[DAQPath folder(i,1).name '\Long_data'])
    else
        disp('No DAQ files found')
    end
    if isempty(dir([IntellivuePath folder(i,1).name '\*.mat']))==0
        copyfile ([IntellivuePath folder(i,1).name '\*.mat'],[IntellivuePath folder(i,1).name '\Long_data'])
    else
        disp('No Intellivue files found')
    end
% Copy split files into real folder overwriting old except for skipped once
    if isempty(dir([camPath folder(i,1).name '\splittet_data\*.avi']))==0
        copyfile ([camPath folder(i,1).name '\splittet_data\*.avi'],[camPath folder(i,1).name])
        delete([camPath folder(i,1).name '\splittet_data\']);
    else
        disp('No splitted videos found')
        delete([camPath folder(i,1).name '\splittet_data\']);        
    end
    if isempty(dir([DAQPath folder(i,1).name '\splittet_data\*.mat']))==0    
        copyfile ([DAQPath folder(i,1).name '\splittet_data\*.mat'],[DAQPath folder(i,1).name])
        delete([DAQPath folder(i,1).name '\splittet_data\']);
    else
        disp('No splittet DAQ files found')
        delete([DAQPath folder(i,1).name '\splittet_data\']);        
    end
    if isempty(dir([IntellivuePath folder(i,1).name '\splittet_data\*.mat']))==0
        copyfile ([IntellivuePath folder(i,1).name '\splittet_data\*.mat'],[IntellivuePath folder(i,1).name])
        delete([IntellivuePath folder(i,1).name '\splittet_data\']);
    else
        disp('No splittet Intellivue files found')
        delete([DAQPath folder(i,1).name '\splittet_data\']);        
    end

end


toc
end