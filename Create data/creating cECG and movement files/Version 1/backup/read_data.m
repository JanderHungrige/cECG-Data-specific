% READ_DATA.M
%
% This file:
% - reads the dataset to be analyzed;
% - is not needed for the online demo on LabView;
%
% This file:
% - is called by launch_processing.m;
% - calls the files main_algorithm.m and offline_analysis.m
%
% aline.serteyn@gmail.com, martijn.schellekens@philips.com
% June 2012




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                      READ DATA                               %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize variables
% --------------------
offsetCAM=0;
offsetPRESS=0;
offsetECG=0;
fileNumberPRESS=0;
fileNumberCAM=0;
firstFramePRESS=0;
if VERBOSE
    counter=0;
end

if VERBOSE && STARTtime~=0 || ENDtime~= 0
    disp(['A data segment of ' num2str(ENDtime-STARTtime) ' minutes has been selected for the processing.'])
end
cd(DataDirectory)

% Read text and binary files
% ---------------------------
if numFiles>=2
    VIDEO=0; %video and pressure mat are only synchronized for the 1st recorded ECG files
    CAMERA=0;
end
for fileNumber=0:numFiles-1
    disp(['Numer of files' num2str(numFiles)]);
    if fileNumber>0
        offsetECG=offsetECG+sizeFile;
    end
    %cd(DataDirectory)
    % Get size of the ECG and RESP files
    % """""""""""""""""""""""""""""""""""
    fid = fopen(['Ref ECG1_0000000',num2str(fileNumber),'.txt'], 'r');
    temp = fread(fid,inf,'*char')';
%     %%%%%% INSERTED BY JAN
                mref = eval(['memmapfile(''Ref ECG1_0000000',num2str(fileNumber),'.bin'',''format'',''uint8'');']);
                dat=mref.data;
                ECG_bin_cell{fileNumber+1}=double(typecast(dat(1:3:end),'int8')')*2^16+double(dat(2:3:end)')*256+double(dat(3:3:end)');  % load all ECG data, +1 as cell has to start with 1 not 0
            % calc total length of analog ECG from data length
                cellsz = cellfun(@sum,cellfun(@numel,ECG_bin_cell,'uni',false),'uni',false);
                sizeFile=sum(cellfun(@sum,cellsz)); clearvars cellsz;
                disp(['J: Size ECG (and RESP) file #' num2str(fileNumber+1) ': ' num2str(sizeFile) ' samples (' num2str(sizeFile/8000/60) ' minutes)']);
    %%%%%% INSERTED BY JAN
    
%   %%% Outcommented by Jan  ************************          
%     a=8; sizeFile=NaN;
%     while isnan(sizeFile)
%         sizeFile=str2double(temp(end-a:end));
%         a=a-1;
%     end
%     timeStartECG=str2double(temp(70:88)); %see text file
%     fclose(fid);
%     
%     if sizeFile==4000
%         if ENDtime>=Inf
%             disp('------>>>ERROR<<<<------: the length of the recording is unknown because the LabView recording software crashed. Please, manually enter the length of the recording (in number of samples) in the Ref ECG1_00000000.txt file!!! (replace 4000 by a larger value)')
%             break
%         else
%             sizeFile=ENDtime*60*8000; %in samples
%         end
%     else
%         if VERBOSE
%             disp(['Size ECG (and RESP) file #' num2str(fileNumber+1) ': ' num2str(sizeFile) ' samples (' num2str(sizeFile/8000/60) ' minutes)']);
%         end
%     end
    % outcommented by Jan ******************************
    
    
    if VIDEO || CAMERA || PRESSMAT
        % Get size of the PRESSURE files
        % """""""""""""""""""""""""""""""
        fid = fopen(['Tekscan_Image_0000000',num2str(fileNumberPRESS),'.txt'], 'r');
        temp = fread(fid,inf,'*char')';
        a=6; sizeFilePRESS=NaN;
        while isnan(sizeFilePRESS)
            sizeFilePRESS=str2double(temp(end-a:end));
            a=a-1;
        end
        timeStartPRESS=str2double(temp(70:88));
        fclose(fid);
        if sizeFilePRESS==1
            %'------>>>ERROR<<<<------: the length of the recording is unknown because the LabView recording software crashed. Please, manually enter the length of the recording (in number of samples) in the Ref ECG1_00000000.txt file!!! (replace 4000 by a larger value)')
            sizeFilePRESS=ENDtime*60*60; %max 60 frames per second?
        end
        % Get size of the CAMERA files
        % """""""""""""""""""""""""""""
        fid = fopen(['uEye_Video_0000000',num2str(fileNumberCAM),'.txt'], 'r');
        temp = fread(fid,inf,'*char')';
        a=6; sizeFileCAM=NaN;
        while isnan(sizeFileCAM)
            sizeFileCAM=str2double(temp(end-a:end));
            a=a-1;
        end
        timeStartCAM=str2double(temp(70:88));
        fclose(fid);
        % Get initial offsets between the different signals
        % """""""""""""""""""""""""""""""""""""""""""""""""
        %expressed in number of ECG samples because the sliding windows for analysis (numberSamples) and the general counter (sampleCounter) counts in term of ECG samples and not seconds nor minutes
        offsetCAM=round((timeStartCAM-timeStartECG)*10^(-7)*8000);
        offsetPRESS=round((timeStartPRESS-timeStartECG)*10^(-7)*8000);
        if VERBOSE
            disp(['The CAMERA recordings and PRESSURE MAT recordings have a delay of, respectively, ' num2str(offsetCAM/8000) ' seconds and ' num2str(offsetPRESS/8000) ' seconds compared to the ECG recordings']);
        end
    else
        offsetCAM=0;
        offsetPRESS=0;
    end
% ******************************************
%                Load data
% ******************************************   
    % """""""""
    if fileNumber==0 %first ECG file
        sampleCounter=max(offsetCAM,offsetPRESS)+STARTtime*60*8000; %define the first ECG sample to be analyzed
        if SAVE
            initOffset=sampleCounter/8000; %intial offset in seconds
        end
    else
        sampleCounter=offsetECG; %assumes that the ECG bin files follow each other perfectly
        if STARTtime*60*8000>offsetECG
            sampleCounter=STARTtime*60*8000; %AS lose synchro CAM/PRESS files
        end
    end
    
    numberNewSamples=WINDOW_SIZE;
    while sampleCounter-offsetECG+numberNewSamples<=sizeFile
        if VIDEO||CAMERA||PRESSMAT
            if firstFramePRESS+numberNewSamples/8000*60 > sizeFilePRESS %max 60 frames per second
                fileNumberPRESS=fileNumberPRESS+1;
                offsetPRESS=sampleCounter;
                if VERBOSE
                    disp(['END OF THE PRESS FILE #' num2str(fileNumberPRESS-1)])
                end
            end
            if (sampleCounter-offsetCAM)/8000*8 >= sizeFileCAM %8 frames per second
                fileNumberCAM=fileNumberCAM+1;
                offsetCAM=sampleCounter;
                if VERBOSE
                    disp(['END OF THE CAM FILE #' num2str(fileNumberCAM-1)])
                end
            end
        end
        %cd(DataDirectory)
% ************** Read reference ECG ******************
        mref = eval(['memmapfile(''Ref ECG1_0000000',num2str(fileNumber),'.bin'',''format'',''uint8'');']);
        xref = mref.data(3*(sampleCounter-offsetECG)+1:3*(sampleCounter-offsetECG)+3*numberNewSamples);
        LengthData = length(xref); % in bytes
        ECG1 = double(typecast(xref(1:3:LengthData),'int8')')*2^16+double(xref(2:3:LengthData)')*256+double(xref(3:3:LengthData)');
        
% ************** Read reference respiration ******************
        %aa=importdata('Resp_01.txt'); aa=aa.data; aa=aa(:,2)';
%         mref = eval(['memmapfile(''Ref Resp_0000000',num2str(fileNumber),'.bin'',''format'',''uint8'');']);
%         xref = mref.data(3*(sampleCounter-offsetECG)+1:3*(sampleCounter-offsetECG)+3*numberNewSamples);
%         LengthData = length(xref); % in bytes
%         referenceRespi = double(typecast(xref(1:3:LengthData),'int8')')*2^16+double(xref(2:3:LengthData)')*256+double(xref(3:3:LengthData)');

% ************** Read capacitive sensors ******************
        for sens=1:8
            m = eval(['memmapfile(''Sensor ',num2str(sens),'_0000000',num2str(fileNumber),'.bin'',''format'',''uint8'');']);
            x = m.data(3*(sampleCounter-offsetECG)+1:3*(sampleCounter-offsetECG)+3*numberNewSamples);
           eval(['SEN',num2str(sens),' = double(typecast(x(1:3:LengthData),''int8'')'')*2^16+double(x(2:3:LengthData)'')*256+double(x(3:3:LengthData)'');']);
        end
        
        %%%%% included by Jan to only save the raw cECG
        Fs= 8000;% smaple frequency of the cECG
        
        information=cell(3,2); %creating an information "txt" file
           information{1,1}='Raw_cECGMAtrix';
           information{1,2}='are all 8 sensors in columns. Each column for one sensor';
           information{2,1}= 'SEN1-8';
           information{2,2}='are the single sensors';
           information{3,1}= 'Fs';
           information{3,3}= 'Sample frequency';
        Raw_cECGMatrix=[SEN1;SEN2;SEN3;SEN4;SEN5;SEN6;SEN7;SEN8];
        if saving_RAW_cECG==1 % added by jan 
           save([dataSetpath{2} 'RAWcECG'],'information','Raw_cECGMatrix','SEN1','SEN2','SEN3','SEN4','SEN5','SEN6','SEN7','SEN8','Fs')
        end
        %%%%included by jan   
        
        if VIDEO || PRESSMAT
% ************** Read pressure mat ******************
            if fileNumberPRESS <= 9
                mpress = eval(['memmapfile(''Tekscan_Image_0000000' num2str(fileNumberPRESS) '.bin'',''format'',{''uint8'' [48 44] ''bdata''});']);
            else
                mpress = eval(['memmapfile(''Tekscan_Image_000000' num2str(fileNumberPRESS) '.bin'',''format'',{''uint8'' [48 44] ''bdata''});']);
            end
            mpresstime =  eval(['memmapfile(''Tekscan_TimeStamp_0000000' num2str(fileNumberPRESS) '.bin'',''format'',''uint8'');']);
            xpresstime = mpresstime.data(1:min(8*(floor(max((sampleCounter-offsetPRESS)/8000*60,0)+numberNewSamples/8000*60)),8*sizeFilePRESS)); %AS: should be read in segments! %64 bits %NB: the sampling frequency of the Tekscan is between 50 and 60 Hz.
            LengthData = length(xpresstime);
            timePRESS =double(typecast(xpresstime(1:8:LengthData),'int8')')*2^56+double(xpresstime(2:8:LengthData)')*2^48+double(xpresstime(3:8:LengthData)')*2^40+double(xpresstime(4:8:LengthData)')*2^32+double(xpresstime(5:8:LengthData)')*2^24+double(xpresstime(6:8:LengthData)')*2^16+double(xpresstime(7:8:LengthData)')*2^8+double(xpresstime(8:8:LengthData)');
            %Find first frame to be displayed (based on timestamps!)
            firstFramePRESS=find(timePRESS*10^(-7)>=timeStartECG*10^(-7)+sampleCounter/8000,1); % comparison in seconds
        end
        if VIDEO || CAMERA
% ************** Camera ******************
            if fileNumberCAM <= 9
                mcam = eval(['memmapfile(''uEye_Video_0000000' num2str(fileNumberCAM) '.bin'',''format'',{''uint8'' [752 480] ''bdata''});']);
            elseif fileNumberCAM <= 99
                mcam = eval(['memmapfile(''uEye_Video_000000' num2str(fileNumberCAM) '.bin'',''format'',{''uint8'' [752 480] ''bdata''});']);
            else
                mcam = eval(['memmapfile(''uEye_Video_00000' num2str(fileNumberCAM) '.bin'',''format'',{''uint8'' [752 480] ''bdata''});']);
            end
            firstFrameCAM=max(floor((sampleCounter-offsetCAM)/1000),0)+1; %read the first frame within the processing windows
            %lastFrameCAM=floor(max((sampleCounter-offsetCAM)/1000,0)+numberNewSamples/1000);
        end
        %cd(CurrentDirectory)
        sampleCounter=sampleCounter+numberNewSamples;
        if sampleCounter >= ENDtime*60*8000
            break
        end
        main_algorithm %<-------------------------------------------------
        numberNewSamples=WINDOW_SIZE;
        if VERBOSE && mod(sampleCounter,8000*60)==0
            counter=counter+1;
            disp(['Already ' num2str(counter) ' minutes of data are processed.'])
        end
    end
    if sampleCounter >= ENDtime*60*8000
        break
    end
end

if VIDEO
    close(vidObj);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                FURTHER PROCESS THE SAVED SIGNALS             %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(CurrentDirectory)
if SAVE==1
    tic
    offline_analysis
    toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%             Remarks about the data reading                   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% - the files 0000, 0001, 0002 do not follow each other exactly... loss of a few ms of data in between!!!! (about 200ms)

% - it is assumes that the CAMERA data are evenly sampled at a frequency of 8 Hz (which is not exactly true)

% - the reading of the pressure mat data could be optimized (speed up and reduce memory usage)

% - max 9 different binary files for ECG and RESP, 99 for pressure mat and 999 for video can be read synchronously using this interface

% - only one video frame is read per analysis window (WINDOW_SIZE) 

% - The size of the file written in .txt do not correspond with the real size of the .bin file ?! (AS:To Be Checked)

% - When the LabView recording software crashes, the end time as well as the size of the recorded signal (number of samples) is not written in the .txt file => you should manualy write a signal length in the .txt file OR write a script to compute the signal length (number of samples) based on the file size

% - Impossible to create a video (read the camera and pressure mat data) for ECG data within the second ECG file! (synchronization problems)



