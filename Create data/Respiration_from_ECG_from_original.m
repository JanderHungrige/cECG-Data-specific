%  In A_synch_... we use ralphs R peak detector to calculate the EDR. But
%  Ralphs R peak detector works ometimes better when the ECG is turned (depending on the Electrode configuration)
% Here we recalculate the EDR with this factor 1 or -1 and save it as the
% new EDR signal. Afterwards it has to be reintegrated into the
% folderstruture again (F_...)
% In this m file we do it automatically for all patients. The m file
% Respiration_from_ECG_from_original_manual
% does it for each individual file which is put in manually.
%  

%jan Werth
%-----------------------------------------------------------------------
function Respiration_from_ECG_from_original(bin,intellivue,saveing,factor,dataSetpath)
    tic
% check if there are all intellivue files
    if exist([dataSetpath{1} '\Synched Data\'], 'file')==7
        load([dataSetpath{1} '\Synched Data\Synchronization_info.mat']);
        if intellivue==1
            if exist('Gapdiff','var')==1 % if gapdiff does not exist, there is no proper intellivue file
            else
                intellivue=0;
                disp('No proper intellivue file exists')
            end
        end
    end
                

    cd(dataSetpath{1})
 

    %% ************** Load Bin Data ****************

    %-------------------------------------------------
    %Analog Data (Ref also intelleview) FS=8kHz
    %-------------------------------------------------

    if bin
    %*********************************************************************
    %ECG *****************************************************************
    %*********************************************************************
        ECG_bin_files=dir('Ref ECG1*.bin'); % find all REf ECG files
        numFiles=size(ECG_bin_files,1);%amount of ECG ref files
        if numFiles ~=0 %if there is at least one file
            for fileNumber=0:numFiles-1
                % Read reference ECG from bin files
                    mref = eval(['memmapfile(''Ref ECG1_0000000',num2str(fileNumber),'.bin'',''format'',''uint8'');']);
                    dat=mref.data;
                    ECG_bin_cell{fileNumber+1}=double(typecast(dat(1:3:end),'int8')')*2^16+double(dat(2:3:end)')*256+double(dat(3:3:end)');  % load all ECG data, +1 as cell has to start with 1 not 0
                % calc total length of analog ECG from data length
                    cellsz = cellfun(@sum,cellfun(@numel,ECG_bin_cell,'uni',false),'uni',false);
                    ECG_bin_length=sum(cellfun(@sum,cellsz)); clearvars cellsz;
                    disp([num2str(numFiles),' ECG bin files where found with a total lenght of ', num2str(ECG_bin_length) ])     
            end 
            %calculate FS from txt file
            cd(dataSetpath{1})
            fid = fopen(['Ref ECG1_00000000.txt'], 'r');
    %         temp = fread(fid,inf,'*char')';
            tline = fgetl(fid); %read in txt line by line
            j=1;
            while ischar(tline)
    %             disp(tline)
                temp{j,1} = strsplit(tline,'='); %separate name and value in cell
                j=j+1;
                tline = fgetl(fid); %go to nextline
            end
            temp{3,1}{1,2}=str2double(temp{3,1}{1,2});
            temp{6,1}{1,2}=str2double(temp{6,1}{1,2});
            if size(temp)>=7 %some txt files do not have the field 
                temp{7,1}{1,2}=str2double(temp{7,1}{1,2});
            end

            if temp{3,1}{1,2}~=[]
    %             timeStartECG=str2double(temp(70:88)); %see text file
                timeStartECG=temp{3,1}{1,2};
            else
                disp('no start time in txt');
                check=1;
            end
            if temp{6,1}{1,2}~=[]
    %             timeEndECG=str2double(temp(155:174));
                timeEndECG=temp{6,1}{1,2};
            else
                disp('no end time in txt')
                check=1;
            end
            if size(temp)>=7 
                if temp{7,1}{1,2}
    %             NRofSamples=str2double(temp(195:end)); 
                    NRofSamples=temp{7,1}{1,2};
                else
                    disp('no samples in txt');
                    check=1;
                end
            end
            if exist('check','var')==0
                t=timeEndECG-timeStartECG; %100nanoseconds
                t=t*10^(-7);%seconds
                FS_bin=NRofSamples/t;
            else
                FS_bin=8000;
                disp('FS_bin set to 8k')
            end
            fclose(fid);
            clearvars fid temp timeStartECG timeEndECG NRofSamples check j


        elseif numFiles ==0
            error('No Ref ECG1 bin files could be found')
        end; clearvars fileNumber mref dat ECG_bin_length

        %Create one long ECG file
        ECG_bin_tot = ECG_bin_cell{1,1};
        if numFiles > 1
            for j=2:numFiles
            ECG_bin_tot=[ECG_bin_tot, ECG_bin_cell{1,j}];
            end
        end; clearvars numFiles ECG_bin_data_cell

        %Downsample to 500Hz (8000/4=2000/4=500)
        ECG_bin_tot_500 = decimate(ECG_bin_tot,4); % if factor greater than 13, decide in smaler peacec for better results
        ECG_bin_tot_500 = decimate(ECG_bin_tot_500,4); ECG_bin_tot_500 = ECG_bin_tot_500'; 
        clearvars ECG_bin_data_tot ECG_bin_cell

        clearvars numFiles fileNumber mref dat Resp_bin_cell Resp_bin_tot
    end
    %% ************** Load Intelleview Data ****************
    %-------------------------------------------------
    % INTELLEVIEW DATA % ECG_Fs=500 Hz ; Resp_Fs=62Hz
    %-------------------------------------------------
    if intellivue 
        cd(dataSetpath{2})
        %--------------ECG-----------
        ECG_intelleview_files=dir('II_*.txt'); % find all REf ECG files
        kb=ECG_intelleview_files.bytes;
        if kb <= 100000 %data smaller than 100kbyte goto II...txt
            ECG_intelleview_files=dir('I_*.txt'); kb=ECG_intelleview_files.bytes;
            if kb <= 100000 %data still smaller than 100kbyte goto III...txt
                ECG_intelleview_files=dir('III_*.txt'); 
            end
        end; clearvars kb


        numFiles=size(ECG_intelleview_files,1);%amount of ECG ref files 
        for fileNumber=1:numFiles
        %read in txt files
            filename=ECG_intelleview_files(fileNumber).name;
            delimiter = '\t';
            startRow = 15; % first few lines are nan. Check if always 15

            formatSpec = '%*s%f%[^\n\r]';
            fileID = fopen(filename,'r');
            textscan(fileID, '%[^\n\r]', startRow-1, 'ReturnOnError', false);
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'ReturnOnError', false);
            fclose(fileID);
            ECG_intelleview_cell{fileNumber} = [dataArray{1:end-1}];
            clearvars filename delimiter startRow formatSpec fileID dataArray ans;
        end
        % calc total length of intelleview ECG
        cellsz = cellfun(@sum,cellfun(@numel,ECG_intelleview_cell,'uni',false),'uni',false);
        ECG_inteleview_length=sum(cellfun(@sum,cellsz)); clearvars cellsz;
        disp([num2str(numFiles),'ECG intelleview files where found with a total length of ',num2str(ECG_inteleview_length),' samples' ]);


       %reate one long ECG file
        ECG_intelleview_tot = ECG_intelleview_cell{1,1};
        if numFiles > 1
            for j=2:numFiles
            ECG_intelleview_tot=[ECG_intelleview_tot; ECG_intelleview_cell{1,j}];
            end
        end; clearvars numFiles 

        disp(['To compare: 500Hz bin length is ' num2str(length(ECG_bin_tot_500))]);
        disp(['To compare: 500Hz int length is ' num2str(length(ECG_intelleview_tot))]);

    end % end of intelleview
    clearvars numFiles fileNumber Resp_intelleview_cell


    %%  ************** Calculate Resiration from ECG (EDR) **************

           %EDR should have same length as ECG
        
        if bin 
           [EDR_bin]=Respiration_from_ECG(factor*ECG_bin_tot_500,500); % calculate respiration from ECG
        end
        if intellivue
           [EDR_intelleview]=Respiration_from_ECG(factor*ECG_intelleview_tot,500); % calculate respiration from ECG
        end

    %% ************** Synch EDR **************

    % first to other ECG, than to camera
    if exist([dataSetpath{1} '\Synched Data\'], 'file')==7
        load([dataSetpath{1} '\Synched Data\Synchronization_info.mat']);
        
        if intellivue==1
            if exist('Gapdiff','var')==1 % pat 4 has intellivue but misses one intellivue file
                if Gapdiff >0 %intelleview gap is earlier
                        %add zeroes or nans at the beginning until gapstart ==
                        pasting(1:Gapdiff,1)=NaN;
                        EDR_intelleview_synched=[pasting;EDR_intelleview];
                    elseif Gapdiff <0 % Intelleview Gap is later
                        %Cut intelleview from the begining as it starts earlier than bin
                        EDR_intelleview_synched=EDR_intelleview(-Gapdiff:end,1);
                    elseif Gapdiff == 0
                        EDR_intelleview_synched=EDR_intelleview;
                end
          

                if offsetCAM_500 >0 %Camera starts later then ECG
                    %add nan or zeroes to ECG, intelleview ECG, Respiration
                    addition(1:offsetCAM_500,1)=NaN;

                    EDR_bin_synched          =[addition'  EDR_bin']';
                    EDR_intelleview_synched  =[addition'  EDR_intelleview_synched']';
                elseif offsetCAM_500 <0 %camera started earlier than ECG
                    %cut the first few samples of ECG, intelleview ECG, Respiration

                    EDR_bin_synched          =EDR_bin(-offsetCAM_500:end,1);
                    EDR_intelleview_synched  =EDR_intelleview_synched(-offsetCAM_500:end,1);
                elseif offsetCAM_500==0
                    EDR_bin_synched=EDR_bin;
                    disp('Cam and ECG are same Please check if data is Ok')
                end
                clearvars offsetCam addition temp
            else %if gapdiff
                 if offsetCAM_500 >0 %Camera starts later then ECG
                    %add nan or zeroes to ECG, intelleview ECG, Respiration
                    addition(1:offsetCAM_500,1)=NaN;

                    EDR_bin_synched          =[addition'  EDR_bin']';
                elseif offsetCAM_500 <0 %camera started earlier than ECG
                    %cut the first few samples of ECG, intelleview ECG, Respiration
                    EDR_bin_synched          =EDR_bin(-offsetCAM_500:end,1);
                elseif offsetCAM_500==0
                    EDR_bin_synched=EDR_bin;
                    disp('Cam and ECG are same Please check if data is Ok')
                 end  
            end%if gapdiff

        else % if not intellivue
             if offsetCAM_500 >0 %Camera starts later then ECG
                %add nan or zeroes to ECG, intelleview ECG, Respiration
                addition(1:offsetCAM_500,1)=NaN;

                EDR_bin_synched          =[addition'  EDR_bin']';
            elseif offsetCAM_500 <0 %camera started earlier than ECG
                %cut the first few samples of ECG, intelleview ECG, Respiration
                EDR_bin_synched          =EDR_bin(-offsetCAM_500:end,1);
            elseif offsetCAM_500==0
                EDR_bin_synched=EDR_bin;
                disp('Cam and ECG are same Please check if data is Ok')
            end
            clearvars offsetCam addition temp    
        end
    else
        disp([dataSetpath{1} ' Synch folder does not exist'])
        return
    end%ifexist synch folder

    %% ************** Save Files **************
    cd(dataSetpath{1});
    savedatapath=([dataSetpath{1} '\Synched Data\']);
    cd([dataSetpath{1} '\Synched Data\'])

    if saveing==1
        FS_ecg = 500;
        save([savedatapath 'Synced_bin_EDR_500b'],'EDR_bin_synched','FS_ecg');
        disp('Synched bin ECG RESP EDR has been saved')

        if intellivue==1 && exist('Gapdiff','var')==1
            save([savedatapath '\Synced_intelleview_EDR_500bHz.mat'],'EDR_intelleview_synched');  
            disp('Synched intelleview (ECG, Resp, EDR) data has been saved')
        end
    end
toc
end