% This function is called if the data does not have any intelleview files.
% That is the case for patient 1-4. Then only the bin files are synched
% with the video and an EDR is created.
% The synched files are saved in the datasetpath in a new folder. 
%
% Input: bin= if binary should be considered (default 1)
%        mov= if movement files should be loaded (defualt 0)
%        cam= if bin videos should be loaded (defualt 0) 
%datasetpath= Where to find the data. Should be declared in sync_data.m
% matlabpath= wher the m files are located. SHould be declared in sync_data.m

%Jan Werth

function synch_only_bin(bin,mov,cam,bin_save,dataSetpath,datacreationPath)


%% ************** Load Bin Data ****************
 
%-------------------------------------------------
%Analog Data (Ref also intelleview) FS=8kHz
%-------------------------------------------------
cd(dataSetpath{1})
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
        temp{7,1}{1,2}=str2double(temp{7,1}{1,2});
        
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
        if temp{7,1}{1,2}
%             NRofSamples=str2double(temp(195:end)); 
            NRofSamples=temp{7,1}{1,2};
        else
            disp('no samples in txt');
            check=1;
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

%*********************************************************************
%RESPIRATION *********************************************************
%*********************************************************************
    Resp_bin_files=dir('Ref Resp_*.bin'); % find all Ref Respiration files
    numFiles=size(Resp_bin_files,1);%amount of Respiration ref files
    if numFiles ~=0 %if there is at least one file / only patients 1-7 have respiration .bin files
        for fileNumber=0:numFiles-1
            % Read reference Respiration from bin files
                mref = eval(['memmapfile(''Ref Resp_0000000',num2str(fileNumber),'.bin'',''format'',''uint8'');']);
                dat=mref.data;
                Resp_bin_cell{fileNumber+1}=double(typecast(dat(1:3:end),'int8')')*2^16+double(dat(2:3:end)')*256+double(dat(3:3:end)');  % load all Respiration data, +1 as cell has to start with 1 not 0
                cellsz = cellfun(@sum,cellfun(@numel,Resp_bin_cell,'uni',false),'uni',false);
                Resp_bin_length=sum(cellfun(@sum,cellsz)); clearvars cellsz;
        end 
        
            %Create one long Respiration file
        Resp_bin_tot = Resp_bin_cell{1,1};
        if numFiles > 1
            for j=2:numFiles
                Resp_bin_tot=[Resp_bin_tot, Resp_bin_cell{1,j}];
            end
        end
        
        %Downsample to 500Hz
        Resp_bin_tot_500 = decimate(Resp_bin_tot,4); % if factor greater than 13, decide in smaler peacec for better results
        Resp_bin_tot_500 = decimate(Resp_bin_tot_500,4); %Resp_bin_tot_500 = Resp_bin_tot_500'; 
        Resp_bin_tot_500 = Resp_bin_tot_500';
        
        disp([num2str(numFiles),' Resp bin files where found with a total lenght of ', num2str(Resp_bin_length) ])
        
    elseif numFiles ==0
        disp('No Ref Resp .bin files could be found');
        cd('C:\Users\310122653\Documents\PhD\cECG Data\Matlab')
        Resp_bin_tot_500(1:length(ECG_bin_tot_500),1)=nan; % create empty Resp file to prevent error
    end
    
%****** Calculate Resiration from ECG (EDR) 
            %EDR should have same length as ECG
     cd(datacreationPath)
    [EDR_bin]=Respiration_from_ECG(ECG_bin_tot_500,500); % calculate respiration from ECG
   
    clearvars numFiles fileNumber mref dat Resp_bin_cell Resp_bin_tot
    
%*********************************************************************    
%MOVEMENT with preassure mat *****************************************
%*********************************************************************
    if mov
        Mov_bin_files=dir('Tekscan_Image_*.bin'); % find all tekscan images files
        numFiles=size(Mov_bin_files,1);%amount of ECG ref files
        if numFiles ~=0 %if there is at least one file 
            for fileNumber=0:numFiles-1
                % Read reference Respiration from bin files
                if numFiles<=9
                    mref = eval(['memmapfile(''Tekscan_Image_0000000' num2str(fileNumber) '.bin'',''format'',{''uint8'' [48 44] ''bdata''});']);
                else
                    mref = eval(['memmapfile(''Tekscan_Image_000000' num2str(fileNumber) '.bin'',''format'',{''uint8'' [48 44] ''bdata''});']);
                end
                    dat=mref.data;
                    Mov_tot_cell{fileNumber+1} =double(typecast(dat(1:8:end),'int8')')*2^56+double(dat(2:8:end)')*2^48+double(dat(3:8:end)')*2^40+double(dat(4:8:end)')*2^32+double(dat(5:8:end)')*2^24+double(dat(6:8:end)')*2^16+double(dat(7:8:end)')*2^8+double(dat(8:8:end)');             
                    cellsz = cellfun(@sum,cellfun(@numel,Mov_tot_Cell,'uni',false),'uni',false);
                    Mov_bin_length=sum(cellfun(@sum,cellsz)); clearvars cellsz;              
                    disp([num2str(numFiles),'Tekscan_Image_ bin files where found with a total lenght of ', num2str(Mov_bin_length) ])     
            end
            %create one long movement file
            Mov_bin_tot = Mov_tot_Cell{1,1};
            if numFiles > 1
                for j=2:numFiles
                Mov_bin_tot=[Resp_bin_tot, Mov_tot_Cell{1,j}];
                end
            end
        
        elseif numFiles ==0
            disp('No Tekscan_Image_ bin files could be found')  ;          
        end; 
    

        
                
% % % % % % %     %downsample to 500Hz
% % % % % % %     ECG_bin_data_tot_500 = decimate(ECG_bin_data_tot,8); % if factor greater than 13, decide in smaler peacec for better results
% % % % % % %     ECG_bin_data_tot_500 = decimate(ECG_bin_data_tot_500,2); ECG_bin_data_tot_500 = ECG_bin_data_tot_500'; 
% % % % % % %  

        clearvars numFiles fileNumber mref dat Mov_tot_cell
     end %end of mov

%*********************************************************************     
%CAMERA **************************************************************
%*********************************************************************
    if cam
        Cam_bin_files=dir('uEye_Video_*.bin'); % find all tekscan images files
        numFiles=size(Cam_bin_files,1);%amount of ECG ref files
        if numFiles ~=0 %if there is at least one file 
            for fileNumber=0:numFiles-1
                if fileNumber <=9
                    %Read reference Respiration from bin files
                    mref = eval(['memmapfile(''uEye_Video_0000000' num2str(fileNumber) '.bin'',''format'',{''uint8'' [752 480] ''bdata''});']);
                elseif  fileNumber <=99
                    mref = eval(['memmapfile(''uEye_Video_000000' num2str(fileNumber) '.bin'',''format'',{''uint8'' [752 480] ''bdata''});']);
                elseif  fileNumber <=999
                    mref = eval(['memmapfile(''uEye_Video_00000' num2str(fileNumber) '.bin'',''format'',{''uint8'' [752 480] ''bdata''});']);
                end
                dat=mref.data;
                Cam_tot_cell{fileNumber+1} =double(typecast(dat(1:8:end),'int8')')*2^56+double(dat(2:8:end)')*2^48+double(dat(3:8:end)')*2^40+double(dat(4:8:end)')*2^32+double(dat(5:8:end)')*2^24+double(dat(6:8:end)')*2^16+double(dat(7:8:end)')*2^8+double(dat(8:8:end)');             
                cellsz = cellfun(@sum,cellfun(@numel,Cam_tot_cell,'uni',false),'uni',false);
                Cam_bin_length=sum(cellfun(@sum,cellsz)); clearvars cellsz;            
            end 
        elseif numFiles ==0
            disp('No uEye_Video_ bin files could be found')   ;        
        end; clearvars numFiles fileNumber     %end fileNumber
             
        %create one long CAM file
        Cam_bin_tot = Cam_tot_Cell{1,1};
        if numFiles > 1
            for j=2:numFiles
            Cam_bin_tot=[Cam_bin_tot, Cam_tot_Cell{1,j}];
            end
        end        
    end %if cam
    
end % end of bin

clearvars temp fid dat dat_with_offet mref Cam_tot_cell


cd(dataSetpath{1})
% Get size of the ECG and RESP files
fid = fopen(['Ref ECG1_00000000.txt'], 'r');
temp = fread(fid,inf,'*char')';
timeStartECG=str2double(temp(70:88)); %see text file
fclose(fid); clearvars fid temp 

% Get size of the CAMERA files
fid = fopen(['uEye_Video_00000000.txt'], 'r');
temp = fread(fid,inf,'*char')';
timeStartCAM=str2double(temp(70:88));
fclose(fid); clearvars fid temp 

%offset expressed in number of seconds and ECG samples 
offsetCam =timeStartCAM-timeStartECG;%100nano seconds
offsetCam_s =offsetCam*10^(-7);%seconds
offsetCAM_8k=round(offsetCam_s*8000); % samples % 10^-7 because of .Net System. DateTime which is accurate on 100 nano seconds. Therefor, second base is 10-7
offsetCAM_500=round(offsetCam_s*500); % samples
disp(['Camera ECG offset is: ',num2str(offsetCAM_500),' samples and ' num2str(offsetCam_s), ' seconds']);
clearvars timeStartECG timeStartECG offsetCAM_8k 

if offsetCam >0 %Camera starts later then ECG
    %add nan or zeroes to ECG, intelleview ECG, Respiration
    addition(1:offsetCAM_500,1)=NaN;
    ECG_bin_synched          =[addition'  ECG_bin_tot_500']';
    Resp_bin_synched         =[addition'  Resp_bin_tot_500']' ;
    EDR_bin_synched          =[addition'  EDR_bin']';
elseif offsetCam <0 %camera started earlier than ECG
    %cut the first few samples of ECG, intelleview ECG, Respiration
    ECG_bin_synched          =ECG_bin_tot_500(-offsetCAM_500:end,1);
    Resp_bin_synched         =Resp_bin_tot_500(-offsetCAM_500:end,1);
    EDR_bin_synched          =EDR_bin(-offsetCAM_500:end,1);
elseif offsetCam==0
    ECG_bin_synched=ECG_bin_tot_500;
    Resp_bin_synched=Resp_bin_tot_500;
    EDR_bin_synched=EDR_bin;
    disp('Cam and ECG are same Please check if data is Ok')
end
clearvars offsetCam addition temp

%%%%%%%%%%%%%%%%%%%%%
    disp(['To compare: 500Hz synched bin length is ' num2str(length(ECG_bin_synched))]);
%%%%%%%%%%%%%%%%%%%%%

%% ************** Save Files **************
cd(dataSetpath{1});
mkdir ('Synched Data') % create synched data folder
savedatapath=([dataSetpath{1} '\Synched Data']);
cd([dataSetpath{1} '\Synched Data'])

save([savedatapath '\Synchronization_info'],'offsetCam_s','offsetCAM_500');


if bin_save
    FS_ecg = 500;
    save([savedatapath '\Synced_bin_ECG_500'],'ECG_bin_synched','FS_ecg')
    save([savedatapath '\Synced_bin_Resp_500'],'Resp_bin_synched','FS_ecg')
    save([savedatapath '\Synced_bin_EDR_500'],'EDR_bin_synched','FS_ecg')

    disp('Synched bin ECG RESP EDR has been saved')
end  


clearvars path dataSetpath savedatapath bin cam numFiles_dataSet j
clearvars t1 t2 t_500 t_62 x1 str acor I
clearvars timeStartCam offsetCAM_500 offsetCam_s offsetCam lag lagDiff mov Fortfahren
clearvars ECG_bin_files 
clearvars Gapstart_bin_idx Gapstart_bin_idx_s

end