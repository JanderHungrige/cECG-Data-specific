function   synch_new_format(bin,intelleviue_save, mov,cam,bin_save,dataSetpath,datacreationPath)
% IN THIS VERSION WE USE THE INFORMATION IN THE TDMS FILE FOR EACH VIDEO TO
% DTERMINE THE STAT OF THE VIDEO AND THEN THE DIFFERENCE BETWEEN THE ECG
% START AND VIDEO START. APPARENTLY THERE IS A MISTAKE IN THIS APROACH. tHE
% SIGNALS ARE NOT SYNCHED.
% WHEN COMPARING THE DATE CREATED TIMESTAMP OF THE FILE AND THE TDMS
% INFORAMTION THERE IS A HUGE DIFFERNEC 
% v2 TRIES TO DEAL WITH THIS WHILE USING THE FILE CREATED TIMESTAMP

% with te cECG dataset they changed the saved format with patients 15
% onwards. Now they use TDMS files to store the Video inforamtion instead
% of txt files. It could be that also the intelliviue infos are saved here
% as they have a struct for respiration. .. not sure

%Adding TDMS read functionpath
addpath('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Create data\Read TDSM');
addpath('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Create data\Read TDSM\tdmsSubfunctions');
addpath('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Create data\Read TDSM\private');


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
        Resp_bin_tot_500=Resp_bin_tot_500';
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
        end
    

        
                
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
        numFiles=size(Cam_bin_files,1);%amount 
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
%% ************** Load Intelleview Data ****************
%-------------------------------------------------
% INTELLEVIEW DATA % ECG_Fs=500 Hz ; Resp_Fs=62Hz
%-------------------------------------------------
   
cd(dataSetpath{2})
%--------------ECG-----------
ECG_intelleview_files=dir('II_*.txt'); % find all REf ECG files
kb=ECG_intelleview_files.bytes;
if kb <= 100000 %data smaller than 100kbyte goto I...txt
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


%-----------------Respiration----------------------

Resp_intelleview_files=dir('Resp_*.txt'); % find all Respiration files
numFiles=size(Resp_intelleview_files,1);%amount of Respiration ref files 
for fileNumber=1:numFiles
    filename=Resp_intelleview_files(fileNumber).name;
    delimiter = '\t';
    startRow = 15;

    formatSpec = '%*s%f%[^\n\r]';
    fileID = fopen(filename,'r');
    textscan(fileID, '%[^\n\r]', startRow-1, 'ReturnOnError', false);
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'ReturnOnError', false);
    fclose(fileID);
    Resp_intelleview_cell{fileNumber} = [dataArray{1:end-1}];
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;
end
cellsz = cellfun(@sum,cellfun(@numel,Resp_intelleview_cell,'uni',false),'uni',false);
Resp_inteleview_length=sum(cellfun(@sum,cellsz)); clearvars cellsz;

   %create one long Respiration file
Resp_intelleview_tot = Resp_intelleview_cell{1,1};
if numFiles > 1
    for j=2:numFiles
    Resp_intelleview_tot=[Resp_intelleview_tot; Resp_intelleview_cell{1,j}];
    end
end    

% calculate respiration also from ECG as the respiration data is often
% clipped. The data is already upsampled
   cd('C:\Users\310122653\Documents\PhD\cECG Data\Matlab')
   [EDR_intelleview]=Respiration_from_ECG(ECG_intelleview_tot,500); % calculate respiration from ECG


clearvars numFiles fileNumber Resp_intelleview_cell

%% ************** Synch ECG and ECG **************
%-------------------------------------------------
% Synchronisation ECG-ECG witch correlation
%-------------------------------------------------
% 
[acor,lag] = xcorr(ECG_bin_tot_500,ECG_intelleview_tot);
[~,I] = max(abs(acor));lagDiff = lag(I);
disp(['xcorr lagg is ' num2str(lag(I)) ' samples and '  num2str(lagDiff/500) ' seconds']);
disp('  ')
if lagDiff<0
    disp('laggdiff is negativ')
    disp(' ')
end

t2=linspace(0,length(ECG_intelleview_tot)/500,length(ECG_intelleview_tot));
% t2=(0:length(ECG_intelleview_tot)-1)/500;
t1=linspace(0,length(ECG_bin_tot_500)/500,length(ECG_bin_tot_500));
% t1 = (0:length(ECG_bin_tot_500)-1)/500;

%If the files have been synched before. Here the old synch is loaded to be
%plottet. If there is asyncholder but no files, e.g. due to abbort, just
%delet the folder
if exist([dataSetpath{1} '\Synched Data'],'file')
synchfolder=dir([dataSetpath{1} '\Synched Data']);
load([dataSetpath{1} '\Synched Data\Synchronization_info.mat'],'Gapstart_bin_idx','Gapstart_intelleview_idx')
else
    Gapstart_bin_idx=0;
    disp('no previouse gap mark. Gap mark set to 0(green *)')
end

figure
    subplot(2,1,1)
    plot(t1,ECG_bin_tot_500);hold on
    plot((abs(lagDiff)/500),ECG_bin_tot_500(abs(lagDiff)),'r*')
    if exist([dataSetpath{1} '\Synched Data'],'file')
        plot((abs(Gapstart_bin_idx)/500),ECG_bin_tot_500(abs(Gapstart_bin_idx)),'g*')
    end
    title('ECG_bin_data_tot_500')
    xlabel('Time (s)')

    subplot(2,1,2)
    plot(t2,ECG_intelleview_tot);hold on 
    if lagDiff <=length(ECG_intelleview_tot)
    plot(abs(lagDiff)/500,ECG_intelleview_tot(abs(lagDiff)),'r*')
    end
    if exist([dataSetpath{1} '\Synched Data'],'file')
        plot((abs(Gapstart_intelleview_idx)/500),ECG_intelleview_tot(abs(Gapstart_intelleview_idx)),'g*')
    end
    title('ECG_intelleview_data_tot_500')
    xlabel('Time (s)')
    
  

str = input('Is this the synchronisation gap? Old=o Yes=Y No=enter or abort (A): ','s');

if str=='y' 
    clearvars Gapstart_bin_idx  Gapstart_intelleview_idx
    Fortfahren=0;
    if lagDiff >=0
    ECG_intelleview_synched = ECG_intelleview_tot(lagDiff:end);
    t2 = (0:length(ECG_intelleview_synched)-1)/500;
    elseif lagDiff<0
    add(1:-lagDiff+1)=nan;
    ECG_intelleview_synched = [add,ECG_intelleview_tot];
    end
    savegappdiff=1;
elseif str=='o'
    Fortfahren=1;

elseif str=='a'
     h=findall(0); %closing figures
     delete(h(2:end));
     delete(findall(0,'Type','figure'))
    return %stop execute
else 
    Fortfahren=1;
    clearvars Gapstart_bin_idx  Gapstart_intelleview_idx    
    str = input('Do you want to switch to manual input or gap search? m/g: ','s');
end
%%                      ECG_ECG manual click
%-------------------------------------------------
% Synchronisation ECG-ECG witch manual click
%------------------------------------------------
  if  Fortfahren==1 && str=='m'
    disp('please find gap in bin ECG')
    
  %BIN  
    figure
        plot(ECG_bin_tot_500);hold on
        title('ECG_bin_data_tot_500');       xlabel('Time (s)')
   % Find Sync point by mouse click   
   disp('----------------------------------------------------------------------------')
   disp('Use normal button clicks to add points.')
   disp('A shift-, right-, or double-click adds a final point and ends the selection.')
   disp('Pressing Return  ends the selection without adding a final point.')
   disp('Pressing Backspace or Delete removes the previously selected point.')
   disp('----------------------------------------------------------------------------')
    str = input('Ready to mark the sync point? press any key : ','s');
    
    [x1,~]=getpts; % get mousclick index and y value
    x1=ceil(x1);    %find index in ECG_bin_data_tot_500 closest to x1
    Gapstart_bin_idx=x1(end,1);
    Gapstart_bin_idx_s=x1(end,1)/500;
 % INTELLEVIEW      
    figure
    plot(ECG_intelleview_tot);hold on 
    title('ECG_intelleview_data_tot_500');    xlabel('Time (s)')
        % Find Sync point by mouse click   
   disp('----------------------------------------------------------------------------')
   disp('Use normal button clicks to add points.')
   disp('A shift-, right-, or double-click adds a final point and ends the selection.')
   disp('Pressing Return  ends the selection without adding a final point.')
   disp('Pressing Backspace or Delete removes the previously selected point.')
   disp('----------------------------------------------------------------------------')
    str = input('Ready to mark the sync point? press any key : ','s');
    
    [x1,~]=getpts; % get mousclick index and y value
    x1=ceil(x1);    %find index in ECG_intelleview_data_tot closest to x1
    Gapstart_intelleview_idx=x1(end,1);
    Gapstart_intelleview_idx_s=x1(end,1)/500;
    str='m';
  end %Fortfahren und m  
%%                      Gap detection
%-------------------------------------------------
% Synchronisation ECG-ECG witch RR-peak 
%------------------------------------------------
  if  Fortfahren==1 && str==('g')

    %*************sync Bin ECG and intelleview ECG ****************************
    % RR distance should be high at gap. compare to xcorr
        % rules to find the correct gap
            %1) Search only in first third of dataset
            %2) The values of the gap (+500 samples in and 500 before end) should not be greater than +-  0.1
            %3) the change during the gap should be 0
            %4) The length of the gap must be maximum 1000 samples !!! doesn`t work
            %5) The gap cannot be at the beginning. Search min 10 RR peaks in for maximum

    %ECG R peak detection
    
    %RALPH
    cd('C:\Users\310122653\Documents\PhD\cECG Data\Matlab\R peak detection')

    hr_max=250;ploting=0;saving=0; % timeline in seconds with 500 Hz FS    
    [ecg_r_peak_idx_bin, ~, ~, ~, ~, RR_bin, ~] = ecg_find_rpeaks(t1, ECG_bin_tot_500, 500, hr_max,ploting,saving); 
    [ecg_r_peak_idx_intelleview, ~, ~, ~, ~, RR_intelleview, ~] = ecg_find_rpeaks(t2, ECG_intelleview_tot, 500, hr_max,ploting,saving);  clearvars hr_max saving ploting
    figure
%         subplot(2,1,1);plot(t1,ECG_bin_tot_500); hold on; plot(ecg_r_peak_idx_bin/500,ECG_bin_tot_500(ecg_r_peak_idx_bin),'r*'); title('Ralph bin');
%         subplot(2,1,2);plot(t2,ECG_intelleview_tot); hold on; plot(ecg_r_peak_idx_intelleview/500,ECG_intelleview_tot(ecg_r_peak_idx_intelleview),'r*'); title('Ralph intelleview');
%     %ROJACKERS
%     [out] = streamingpeakdetection(ECG_bin_tot_500', 500, [36 230], 1); %(inputData, Fs, HRLimits, PLOT, Fc, blockSize)
%     ecg_r_peak_idx_bin=out.peakPositionArray;RR_bin=diff(out.peakPositionArray);
%     [out] = streamingpeakdetection(ECG_intelleview_tot', 500, [36 230], 1); %(inputData, Fs, HRLimits, PLOT, Fc, blockSize)
%     ecg_r_peak_idx_intelleview=out.peakPositionArray;RR_intelleviewdiff(out.peakPositionArray);
%     
    %DETERINE GAP IN BIN
    figure
    subplot (2,1,2)
     plot(RR_bin);hold on;xlabel('RR interval lengths for bin')
    subplot (2,1,1)
     plot(ECG_bin_tot_500);hold on
    %  RR_bin=RR_bin(1,1:length(RR_bin)/3);  %1)
    for j=1:length(RR_bin)
        [val,idx]=max(RR_bin); %1) %find the biggest RR interval. This is the gap in the ECG for synchronisation
        if  idx <= floor(length(RR_bin)/3) ...  ;                                                                           %1)
            && ecg_r_peak_idx_bin(idx)-ecg_r_peak_idx_bin(idx-1)>=1000 ...                                                  %4)
            && max(ECG_bin_tot_500(ecg_r_peak_idx_bin(idx-1)+300:(ecg_r_peak_idx_bin(idx)-300))) <= 0.1 ...            %2)
            && max(ECG_bin_tot_500(ecg_r_peak_idx_bin(idx-1)+300:(ecg_r_peak_idx_bin(idx)-300))) >= -0.1 ...           %2)          
            && round(mean(diff(ECG_bin_tot_500(ecg_r_peak_idx_bin(idx-1)+300:(ecg_r_peak_idx_bin(idx)-300)))))==0 ...  %3)
            && idx>10                                                                                                       %5)
                Gapstart_bin_idx=ecg_r_peak_idx_bin(idx-1); %index of the Gap start
                  plot(idx,val,'g*')
                  str = input('Is this the synchronisation gap? y/enter: ','s');
                if str=='y' 
                    break % if statement is reached break the for loop
                else
                     RR_bin(idx)=[];                
                end
        else
            RR_bin(idx)=[]; %delet maximum RR to reach second highest maximum
            for i=1:1 %plotting
                subplot (2,1,1)
                    plot(ecg_r_peak_idx_bin(idx-1),ECG_bin_tot_500(ecg_r_peak_idx_bin(idx-1)),'r*')
                    title('bin ECG')

                subplot (2,1,2)
                    plot(idx,val,'r*')
                    title(['1)' num2str(idx <= floor(length(RR_bin)/3)), ...
                        ' +2)' num2str(max(ECG_bin_tot_500(ecg_r_peak_idx_bin(idx-1)+300:(ecg_r_peak_idx_bin(idx)-300))) <= 0.1), ... 
                        ' -2)' num2str(max(ECG_bin_tot_500(ecg_r_peak_idx_bin(idx-1)+300:(ecg_r_peak_idx_bin(idx)-300))) >= -0.1), ...
                        ' 3)' num2str(max(ECG_bin_tot_500(ecg_r_peak_idx_bin(idx-1)+300:(ecg_r_peak_idx_bin(idx)-300))) >= -0.1), ...
                        ' 4)' num2str(ecg_r_peak_idx_bin(idx)-ecg_r_peak_idx_bin(idx-1)>=1000), ...
                        ' 5)' num2str(idx>10) ...                
                        ]);
            end

            str = input('Is this the synchronisation gap? y/enter: ','s');
            if str=='y'
                disp('This gap is choosen for synchronisation')
                Gapstart_bin_idx=ecg_r_peak_idx_bin(idx-1); %index of the Gap start
                break
            end 
        end % if idx
    end % for

    %DETERMINE GAP IN INTELLEVIEW
    figure
    subplot (2,1,2)
     plot(RR_intelleview);hold on;title('RR interval lengths for intelleview')
     subplot(2,1,1)
     plot(ECG_intelleview_tot);hold on; title('intelleview ECG cut at bin ECG length')
      RR_intelleview=RR_intelleview(1,1:length(RR_bin)/3);  %1)
    for j=1:length(RR_intelleview)
        [val,idx_i]=max(RR_intelleview); %5) %find the biggest RR interval. This is the gap in the ECG for synchronisation
        if  idx_i <= floor(length(RR_intelleview)/3)  ...                                                                                                      %1)
            && ecg_r_peak_idx_intelleview(idx)-ecg_r_peak_idx_intelleview(idx-1)>=1000 ... %4)           
            && max(ECG_intelleview_tot(ecg_r_peak_idx_intelleview(idx_i-1)+300:ecg_r_peak_idx_intelleview(idx_i)-300)) <= 0.1 ...                  %2)
            && max(ECG_intelleview_tot(ecg_r_peak_idx_intelleview(idx_i-1)+300:ecg_r_peak_idx_intelleview(idx_i)-300)) >= -0.1 ...                 %2)
            && round(mean(diff(ECG_intelleview_tot(ecg_r_peak_idx_intelleview(idx_i-1)+300:(ecg_r_peak_idx_intelleview(idx_i-1)-300)))))==0 ...    %3)
            && idx_i>10                                                                                                                                   %5)
                Gapstart_intelleview_idx=ecg_r_peak_idx_intelleview(idx_i-1); %index of the Gap start
                plot(find(RR_intelleview==max(RR_intelleview)),RR_intelleview(find(RR_intelleview==max(RR_intelleview))),'g*')
                  str = input('Is this the synchronisation gap? y/enter: ','s');
                if str=='y' 
                    break % if statement is reached break the for loop
                else
                     RR_intelleview(idx_i)=[];                
                end
        else
            RR_intelleview(idx_i)=[]; %delet maximum RR to reach second highest maximum           
            for i=1:1 % plotting 
                subplot(2,1,1)
                plot(ecg_r_peak_idx_intelleview(idx_i-1),ECG_intelleview_tot(ecg_r_peak_idx_intelleview(idx_i-1)),'r*')

                subplot(2,1,2)
                plot(idx,val,'r*')
                    title(['1)' num2str(idx_i <= floor(length(RR_bin)/3)), ...
                        ' +2)' num2str(max(ECG_intelleview_tot(ecg_r_peak_idx_intelleview(idx_i-1)+300:(ecg_r_peak_idx_intelleview(idx_i)-300))) <= 0.1), ... 
                        ' -2)' num2str(max(ECG_intelleview_tot(ecg_r_peak_idx_intelleview(idx_i-1)+300:(ecg_r_peak_idx_intelleview(idx_i)-300))) >= -0.1), ...
                        ' 3)' num2str(max(ECG_intelleview_tot(ecg_r_peak_idx_intelleview(idx_i-1)+300:(ecg_r_peak_idx_intelleview(idx_i)-300))) >= -0.1), ...
                        ' 4)' num2str(ecg_r_peak_idx_intelleview(idx_i)-ecg_r_peak_idx_intelleview(idx_i-1)>=1000), ...
                        ' 5)' num2str(idx_i>10) ...                
                        ]);
            end
            str = input('Is this the synchronisation gap? y/enter: ','s');
            if str=='y'
                disp('This gap is choosen for synchronisation')
                Gapstart_intelleview_idx=ecg_r_peak_idx_intelleview(idx_i-1); %index of the Gap start
                break
            end
        end %idx_i
    end %for
    str='g'; 
 end % Fortfahren und g

%% ECG-ECG SYNCHRONISATION 
% The data is extended or cut always at the intelleview ECG
    if Fortfahren==1
        Gapdiff=Gapstart_bin_idx-Gapstart_intelleview_idx; % getting data either from manual or gap find loop
        if Gapdiff >0 %intelleview gap is earlier
            %add zeroes or nans at the beginning until gapstart ==
            pasting(1:Gapdiff,1)=NaN;
            ECG_intelleview_synched=[pasting;ECG_intelleview_tot];
        elseif Gapdiff <0 % Intelleview Gap is later
            %Cut intelleview form the begining as it starts earlier than bin
            ECG_intelleview_synched=ECG_intelleview_tot(-Gapdiff:end,1);
        elseif Gapdiff == 0
            ECG_intelleview_synched=ECG_intelleview_tot;
        end
        disp(['The difference between gap_bin and gab_intelle start is ', num2str(Gapdiff), ' samples and ', num2str(Gapdiff/500), ' seconds' ])

    for j=1:1 % Plotting; for loop just to be able to fold the loop together to save space
        figure
       subplot (3,1,1)
        plot(ECG_bin_tot_500);hold on
        plot(Gapstart_bin_idx, ECG_bin_tot_500(Gapstart_bin_idx,1),'r*')
        title('bin ECG')
        xlim([0 length(ECG_bin_tot_500)])
        
       subplot(3,1,2)
        plot(ECG_intelleview_synched);hold on
        title('sync Intelleview ECG ')
        xlim([0 length(ECG_bin_tot_500)])
        
       subplot (3,1,3)
        plot(ECG_intelleview_tot);hold on
        plot(Gapstart_intelleview_idx,ECG_intelleview_tot(Gapstart_intelleview_idx),'r*')
        title('intelleview ECG')
        xlim([0 length(ECG_bin_tot_500)])
        
% 
        if str =='g'
                figure
                plot(RR_intelleview);hold on
                plot(find(RR_intelleview==max(RR_intelleview)),RR_intelleview(find(RR_intelleview==max(RR_intelleview))),'r*')
                title('RR interval lengths for intelleview')

                figure
                plot(RR_bin);hold on
                plot(find(RR_bin==max(RR_bin)),RR_bin(find(RR_bin==max(RR_bin))),'r*')
                title('RR bin lengths for intelleview')
        end
    end 
    
   end % if fortfahre ==1 for synchronisation



clearvars RR_bin idx idx_i RR_intelleview val
%% ************** Synch ECG and Resp **************
%***************sync intelleview ECG and intelleview Respiration ***********
% respiration data FS =62.5 Hz
FS_resp=Resp_inteleview_length/(ECG_inteleview_length/500);
FS_factor=500/FS_resp;
t_62=linspace(0,ceil(Resp_inteleview_length/FS_resp), Resp_inteleview_length)'; % timeline in seconds with 62.5 Hz FS
t_500=linspace(0,ceil(Resp_inteleview_length/FS_resp), Resp_inteleview_length*FS_factor)'; % timeline in seconds with 500Hz FS

Resp_intelleview_tot_500=interp1(t_62,Resp_intelleview_tot,t_500,'pchip');%Shape-preserving piecewise cubic interpolation
disp('')
disp(['To compare: 500Hz intelleview ECG length is  ' num2str(length(ECG_intelleview_tot))]);
disp(['To compare: 500Hz intelleview RESP length is ' num2str(length(Resp_intelleview_tot_500))]);   
Resp_ECG_diff=(length(Resp_intelleview_tot_500)-length(ECG_intelleview_tot));  
disp(['To compare: 500Hz intelleview RESP ECG differnce in length is ' num2str(Resp_ECG_diff)]);

%due to interpolation the Respiation signal can be a bit shorter or longer
%than the ECG. Mostly by one value.
if Resp_ECG_diff < 0
    temp(1:abs(Resp_ECG_diff))=nan;
    Resp_intelleview_tot_500=[Resp_intelleview_tot_500' temp' ]' ;
elseif Resp_ECG_diff > 0
    temp(1:abs(Resp_ECG_diff))=nan;
    Resp_intelleview_tot_500=Resp_intelleview_tot_500(1:end-abs(Resp_ECG_diff));
end
clearvars temp

if abs(length(Resp_intelleview_tot_500)-length(ECG_intelleview_tot))<=5 % if length is approximately the same (can be diff due to interp1)
    if Gapdiff >0 %intelleview gap is earlier
        %add zeroes or nans at the beginning until gapstart ==
        pasting(1:Gapdiff,1)=NaN;
        Resp_intelleview_synched=[pasting;Resp_intelleview_tot_500];
        EDR_intelleview_synched=[pasting;EDR_intelleview];
    elseif Gapdiff <0 % Intelleview Gap is later
        %Cut intelleview from the begining as it starts earlier than bin
        Resp_intelleview_synched=Resp_intelleview_tot_500(-Gapdiff:end,1);
        EDR_intelleview_synched=EDR_intelleview(-Gapdiff:end,1);
    elseif Gapdiff == 0
        Resp_intelleview_synched=Resp_intelleview_tot_500;
        EDR_intelleview_synched=EDR_intelleview;
    end
else
    disp('--------------------------------------------------------------------------------')
    disp('!!Intelleview Respiration and Intelleview ECG are majorly differnt in  length !!')
    disp('--------------------------------------------------------------------------------')

end
    
% for j=1:1 % Plotting; for loop just to be able to fold the loop together to save space
%      figure   
%        subplot(4,1,1)
%         plot(ECG_intelleview_synched);hold on
%         title('sync Intelleview ECG ')
%         xlim([0 length(ECG_bin_data_tot_500)])
%         
%       subplot (4,1,2)
%         plot(Resp_intelleview_synched,1);hold on
%         title('intelleview Resp synched')
%         xlim([0 length(ECG_bin_data_tot_500)])
%         
%        subplot (4,1,3)
%         plot(ECG_intelleview_data_tot);hold on
%         plot(Gapstart_intelleview_idx,ECG_intelleview_data_tot(Gapstart_intelleview_idx),'r*')
%         title('intelleview ECG')
%         xlim([0 length(ECG_bin_data_tot_500)])     
% 
%        subplot (4,1,4)
%         plot(Resp_intelleview_data_tot_500);hold on
%         plot(Gapstart_intelleview_idx,Resp_intelleview_data_tot_500(Gapstart_intelleview_idx),'r*')
%         title('total intelleview Resp')
%    
% end 
clearvars Resp_FS FS_factor 

%% ************** Synch ECG(bin), ECG(intelleview), Resp(intelleview) and Video **************
%****************Sync Video and .bin ECG **********************************
cd(dataSetpath{1})
% Get size of the ECG and RESP files
fid = fopen(['Ref ECG1_00000000.txt'], 'r');
temp = fread(fid,inf,'*char')';
timeStartECG=str2double(temp(70:88)); %str2num(['int64(' temp(70:88) ')']);%see text file
fclose(fid); clearvars fid temp 

% Get size of the CAMERA files trhough the TDMS file

% change the name of the DAQ folder to see if it fits to a video name
proxname=split(dataSetpath{1},'_');proxname=proxname{end,1};
proxnameminus=str2num(proxname)-100; proxnameminus=num2str(proxnameminus); % As the video and DAQ can change in time create also +- one minute
proxnameplus=str2num(proxname)+100; proxnameplus=num2str(proxnameplus);

proxname=[proxname(1,1:2) '.' proxname(1,3:4) '.']; % Make it look the same as in the video folder
videonames=dir([dataSetpath{3} '*' proxname '*Acquisition.tdms']);
if isempty(videonames)==1
    proxname=[proxnameminus(1,1:2) '.' proxnameminus(1,3:4) '.']; % Make it look the same as in the video folder
    videonames=dir([dataSetpath{3} '*' proxname '*Acquisition.tdms']);
end
if isempty(videonames)==1
    proxname=[proxnameplus(1,1:2) '.' proxnameplus(1,3:4) '.']; % Make it look the same as in the video folder
    videonames=dir([dataSetpath{3} '*' proxname '*Acquisition.tdms']);
    videonIDXames=dir([dataSetpath{3} '*' proxname '*Acquisition.tdms_index']);
end
  
% loading camera start time datenum format
%START ORIGINAL
% [finalOutput,metaStruct] = TDMS_readTDMSFile([dataSetpath{3} videonames.name] );% load the tdms file with the same name as the DAQ folder
% Index = find(contains(finalOutput.groupNames,'Instr_uEyeAcquisition')); % get the indx of the uEYE timestamp from the struct
% camtimeIDX=finalOutput.chanIndices{1,Index}; % Get the field index of the Cam time 
% timeStartCam=(finalOutput.data{1,camtimeIDX}(1,1)); %(1,1) is the first time, therfore 
% % Adding the 7 hours to cam time
% timeStartCam = (datetime(timeStartCam, 'ConvertFrom', 'datenum')); % create datetime for easier manipulation https://stackoverflow.com/questions/44158808/time-stamp-conversion-to-seconds#
% format long
% timeStartCam.Hour = timeStartCam.Hour + 7; % correct for the offset of 7 hours. ?timezone? format long
% timeStartCam = datenum(timeStartCam); % create datetime for easier manipulation https://stackoverflow.com/questions/44158808/time-stamp-conversion-to-seconds#
% % Getting the two time systems on the same setting
% % explanation:
% %https://stackoverflow.com/questions/44210768/comparing-datenum-and-18-digit-timestamps
% %check ECGtime
% cehktimeECG=(datetime(timeStartECG, 'ConvertFrom', 'datenum'));
% 
% timeStartECG = int64(timeStartECG);   % .NET datetime number. Change to int64 to not loose precicion
% timeStartCam = timeStartCam;          % MATLAB datetime number
% ticks1_ECG = System.DateTime(timeStartECG).AddYears(1).AddDays(2).Ticks; % add the difference between time systems starting points in Ticks (100NS) to .NET 
% % %original
% ticks2_CAM = int64(timeStartCam*1e9*24*36);%+7*3600; % convert to Ticks (100ns)
% % %orignial
% ticks2_CAM2 = int64(timeStartCam)*1e7*24*3600;%+7*3600; % convert to Ticks (100ns) / protion of the day*24=hours *3600=seconds *1e7=100ns
% 
% dTicks = ticks2_CAM-ticks1_ECG;       % Difference in ticks (100 nanoseconds)
% offsetCam_s = double(dTicks)*(1e-7);  % Difference in seconds. 10^-7 because of .Net System time which is accurate on 100 nano seconds. Therefore, second based is 10-7
%END ORIGINAL

[finalOutput,metaStruct] = TDMS_readTDMSFile([dataSetpath{3} videonames.name] );% load the tdms file with the same name as the DAQ folder
Index = find(contains(finalOutput.groupNames,'Instr_uEyeAcquisition')); % get the indx of the uEYE timestamp from the struct
camtimeIDX=finalOutput.chanIndices{1,Index}; % Get the field index of the Cam time 
timeStartCam=(finalOutput.data{1,camtimeIDX}(1,1)); %(1,1) is the first time, therfore 
% Adding the 7 hours to cam time
timeStartCam = (datetime(timeStartCam, 'ConvertFrom', 'datenum')); % create datetime for easier manipulation https://stackoverflow.com/questions/44158808/time-stamp-conversion-to-seconds#
format long
timeStartCam.Hour = timeStartCam.Hour + 7; % correct for the offset of 7 hours. ?timezone? format long
timeStartCam = datenum(timeStartCam); % .NET datetime number.% create datetime for easier manipulation https://stackoverflow.com/questions/44158808/time-stamp-conversion-to-seconds#
% Getting the two time systems on the same setting
% explanation:
%https://stackoverflow.com/questions/44210768/comparing-datenum-and-18-digit-timestamps
%check ECGtime
% cehktimeECG=(datetime(timeStartECG, 'ConvertFrom', 'datenum'));
ticks1_ECG = (System.DateTime(timeStartECG).AddYears(1).AddDays(2).Ticks)*1e-4; %1ms /// add the difference between time systems starting points in Ticks (100NS) to .NET 
% %original
% ticks2_CAM = (timeStartCam*1e5*24*36);%+7*3600; % convert to Ticks (100ns)
% %orignial
ticks2_CAM = (timeStartCam*1e3*24*3600);%+7*3600; % convert to Ticks (100ns) / protion of the day*24=hours *3600=seconds *1e7=100ns

dTicks = ticks2_CAM-ticks1_ECG;       % Difference in ticks (100 nanoseconds)
offsetCam_s = double(dTicks)*(1e-3);  % Difference in seconds. 10^-7 because of .Net System time which is accurate on 100 nano seconds. Therefore, second based is 10-7

%offset expressed in number of ECG samples 
offsetCAM_8k=round(offsetCam_s*8000); % samples % 
offsetCAM_500=round(offsetCam_s*500); % samples
disp(['Camera ECG offset is: ',num2str(offsetCAM_500),' samples and ' num2str(offsetCam_s), ' seconds']);
clearvars timeStartECG timeStartECG offsetCAM_8k 

if offsetCam_s >0 %Camera starts later then ECG
    %add nan or zeroes to ECG, intelleview ECG, Respiration
    addition(1:offsetCAM_500,1)=NaN;
    ECG_bin_synched          =[addition'  ECG_bin_tot_500']';
    Resp_bin_synched         =[addition'  Resp_bin_tot_500']' ;
    EDR_bin_synched          =[addition'  EDR_bin']';
    Resp_intelleview_synched =[addition'  Resp_intelleview_synched']';
    EDR_intelleview_synched  =[addition'  EDR_intelleview_synched']';
    ECG_intelleview_synched  =[addition'  ECG_intelleview_synched']';
elseif offsetCam_s <0 %camera started earlier than ECG
    %cut the first few samples of ECG, intelleview ECG, Respiration
    ECG_bin_synched          =ECG_bin_tot_500(abs(offsetCAM_500):end,1);
    Resp_bin_synched         =Resp_bin_tot_500(abs(offsetCAM_500):end,1);
    EDR_bin_synched          =EDR_bin(abs(offsetCAM_500):end,1);
    Resp_intelleview_synched =Resp_intelleview_synched(abs(offsetCAM_500):end,1);
    EDR_intelleview_synched  =EDR_intelleview_synched(abs(offsetCAM_500):end,1);
    ECG_intelleview_synched  =ECG_intelleview_synched(abs(offsetCAM_500):end,1);
elseif offsetCam_s==0
    ECG_bin_synched=ECG_bin_tot_500;
    Resp_bin_synched=Resp_bin_tot_500;
    EDR_bin_synched=EDR_bin;
    disp('Cam and ECG are same Please check if data is Ok')
end

%ORIGINAL
% if offsetCam_s >0 %Camera starts later then ECG
%     %add nan or zeroes to ECG, intelleview ECG, Respiration
%     addition(1:offsetCAM_500,1)=NaN;
%     ECG_bin_synched          =[addition'  ECG_bin_tot_500']';
%     Resp_bin_synched         =[addition'  Resp_bin_tot_500']' ;
%     EDR_bin_synched          =[addition'  EDR_bin']';
%     Resp_intelleview_synched =[addition'  Resp_intelleview_synched']';
%     EDR_intelleview_synched  =[addition'  EDR_intelleview_synched']';
%     ECG_intelleview_synched  =[addition'  ECG_intelleview_synched']';
% elseif offsetCam_s <0 %camera started earlier than ECG
%     %cut the first few samples of ECG, intelleview ECG, Respiration
%     ECG_bin_synched          =ECG_bin_tot_500(-offsetCAM_500:end,1);
%     Resp_bin_synched         =Resp_bin_tot_500(-offsetCAM_500:end,1);
%     EDR_bin_synched          =EDR_bin(-offsetCAM_500:end,1);
%     Resp_intelleview_synched =Resp_intelleview_synched(-offsetCAM_500:end,1);
%     EDR_intelleview_synched  =EDR_intelleview_synched(-offsetCAM_500:end,1);
%     ECG_intelleview_synched  =ECG_intelleview_synched(-offsetCAM_500:end,1);
% elseif offsetCam_s==0
%     ECG_bin_synched=ECG_bin_tot_500;
%     Resp_bin_synched=Resp_bin_tot_500;
%     EDR_bin_synched=EDR_bin;
%     disp('Cam and ECG are same Please check if data is Ok')
% end
clearvars offsetCam ticks1_ECG ticks2_ CAM addition temp dticks

%%%%%%%%%%%%%%%%%%%%%
    disp(['To compare: 500Hz synched bin length is ' num2str(length(ECG_bin_synched))]);
    disp(['To compare: 500Hz synched int length is ' num2str(length(ECG_intelleview_synched))]);
%%%%%%%%%%%%%%%%%%%%%


%%  *********** Find the ECG signal which match in length with Video and create same length

% Video_bin_files=dir('uEye_TimeStamp_*.txt'); % find all REf ECG files
% numFiles=size(Video_bin_files,1);%amount of ECG ref files
% % Get size of the CAMERA files. Load startTime form first file and endTime from last file
% lastVideofile=Video_bin_files(numFiles,1).name;
% fid = fopen(lastVideofile, 'r');
% temp = fread(fid,inf,'*char')';
% timeEndCAM=str2double(temp(155:174));
% Videolength=timeEndCAM-timeStartCAM;%100nano seconds
% Videolength_s=Videolength*10^(-7);%seconds
% Videolength_500=round(Videolength_s*500); % sample
% 
% v = VideoReader([lastVideofile '.avi'])
% 
% lengthdiff_between_bin_Cam=abs(length(ECG_bin_synched)-Videolength_500);
% lengthdiff_between_intelleview_Cam=abs(length(ECG_intelleview_synched)-Videolength_500);
% clearvars temp
% 
% if lengthdiff_between_bin_Cam < lengthdiff_between_intelleview_Cam
%         disp('the bin ECG is clostest in length to the Video')
%         This_ECG_and_Video_match=('bin ECG');
%         if length(ECG_bin_synched)>Videolength_500 % if the bin ECG is longer than the Video cut it at the end
%             ECG_bin_synched=ECG_bin_synched(1:Videolength_500,:);
%             EDR_bin_synched=EDR_bin_synched(1:Videolength_500,:);
%             Resp_bin_synched=Resp_bin_synched(1:Videolength_500,:);
%         elseif length(ECG_bin_synched)<Videolength_500 % if bin ECG is shortr than Video, add zeroes at the end
%             temp(1:lengthdiff_between_bin_Cam,1)=nan;
%             ECG_bin_synched=[ECG_bin_synched' temp']';
%             EDR_bin_synched=[EDR_bin_synched' temp']';
%             Resp_bin_synched=[Resp_bin_synched' temp']';clearvars temp
%         elseif length(ECG_bin_synched)==Videolength_500
%             disp('bin ECG and video are the same length')
%         end
%         
% elseif lengthdiff_between_bin_Cam > lengthdiff_between_intelleview_Cam
%         disp('the intelleview ECG is clostest in length to the Video')
%         This_ECG_and_Video_match=('intelleview ECG');        
%         if length(ECG_intelleview_synched)>Videolength_500 % if the intelleview ECG is longer than the Video cut it at the end
%             ECG_intelleview_synched=ECG_intelleview_synched(1:Videolength_500,:); %same with respiration as the interpolated resp is the same length as EcG
%             Resp_intelleview_synched=Resp_intelleview_synched(1:Videolength_500,:);   
%             EDR_intelleview_synched=EDR_intelleview_synched(1:Videolength_500,:);
%         elseif length(ECG_intelleview_synched)<Videolength_500 % if intelleview ECGis shorter than Video, add zeroes at the end
%             temp(1:lengthdiff_between_intelleview_Cam,1)=nan;
%             ECG_intelleview_synched=[ECG_intelleview_synched' temp']';
%             Resp_intelleview_synched=[Resp_intelleview_synched' temp']';
%             EDR_intelleview_synched=[EDR_intelleview_synched' temp']'; 
%             clearvars temp  
%         elseif length(ECG_intelleview_synched)==Videolength_500
%             disp('intelleview ECG, respiration and video are the same length')
%         end
%         
% elseif lengthdiff_between_bin_Cam == lengthdiff_between_intelleview_Cam
%         disp('both ECGs are clostest in length to the Video')
%         This_ECG_and_Video_match=('intelleview and bin ECG');        
%         if length(ECG_intelleview_synched)>Videolength_500 % if the bin ECG is longer than the Video cut it at the end
%             ECG_intelleview_synched=ECG_intelleview_synched(1:Videolength_500,:);
%             Resp_intelleview_synched=Resp_intelleview_synched(1:Videolength_500,:);   
%             EDR_intelleview_synched=EDR_intelleview_synched(1:Videolength_500,:);
%             ECG_bin_synched=ECG_bin_synched(1:Videolength_500,:);
%             EDR_bin_synched=EDR_bin_synched(1:Videolength_500,:);
%             Resp_bin_synched=Resp_bin_synched(1:Videolength_500,:);            
%         elseif length(ECG_intelleview_synched)<Videolength_500 % if bin ECGis shortr than Video, add zeroes at the end
%             temp(1:lengthdiff_between_intelleview_Cam,1)=nan;
%             ECG_intelleview_synched=[ECG_intelleview_synched' temp']';
%             Resp_intelleview_synched=[Resp_intelleview_synched' temp']';
%             EDR_intelleview_synched=[EDR_intelleview_synched' temp']';  clearvars temp           
%             temp(1:lengthdiff_between_bin_Cam,1)=nan;
%             ECG_bin_synched=[ECG_bin_synched' temp'];
%             EDR_bin_synched=[EDR_bin_synched' temp'];
%             Resp_bin_synched=[Resp_bin_synched' temp'];clearvars temp            
%         elseif length(ECG_intelleview_synched)==Videolength_500
%             disp('intelleview and bin ECGs and the video are the same length')
%         end
%         
% end
% fclose(fid); clearvars fid temp 

%% ************** Save Files **************
cd(dataSetpath{1});
mkdir ('Synched Data') % create synched data folder
savedatapath=([dataSetpath{1} '\Synched Data']);
cd([dataSetpath{1} '\Synched Data'])

if exist('Gapdiff','var')
    save([savedatapath '\Synchronization_info'],'Gapdiff','Gapstart_bin_idx','Gapstart_intelleview_idx','offsetCam_s','offsetCAM_500');
end

if bin_save
    FS_ecg = 500;
    save([savedatapath '\Synced_bin_ECG_500'],'ECG_bin_synched','FS_ecg')
    save([savedatapath '\Synced_bin_Resp_500'],'Resp_bin_synched','FS_ecg')
    save([savedatapath '\Synced_bin_EDR_500'],'EDR_bin_synched','FS_ecg')

    disp('Synched bin ECG RESP EDR has been saved')
end
if intelleviue_save
    save([savedatapath '\Synced_intelleview_ECG_500Hz'],'ECG_intelleview_synched','FS_ecg');
    save([savedatapath '\Synced_intelleview_Resp_500Hz'],'Resp_intelleview_synched')  ;  
    save([savedatapath '\Synced_intelleview_EDR_500Hz'],'EDR_intelleview_synched')  ;  
    
    disp('Synched intelleview (ECG, Resp, EDR) data has been saved')
end
%     save([savedatapath '\This_ECG_and_Video_match'],'This_ECG_and_Video_match');
    


clearvars path dataSetpath savedatapath bin intelleview bin_save intelleview_save cam numFiles_dataSet j
clearvars t1 t2 t_500 t_62 x1 str acor I
clearvars timeStartCam offsetCAM_500 offsetCam_s offsetCam lag lagDiff mov Fortfahren
clearvars ECG_bin_files ECG_intelleview_data_cell
clearvars Gapstart_bin_idx Gapstart_bin_idx_s Gapstart_intelleview_idx Gapstart_intelleview_idx_s
%% close figures
input('press any key to close all Figures: ','s');

 h=findall(0);
 delete(h(2:end));

delete(findall(0,'Type','figure'))
clearvars h

toc