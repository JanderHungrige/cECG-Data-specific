% FInd VIdeo ECG shift with annnotations
% same as the read_annotation_csv.m but with plot and loading of ECG files 
%Reading in the annotation file for the cECG dataset
clc
clear
tic

pat=14;
epoch=1; %[30] [1] 30 seconds or 1 secondepoch annotations

addpath('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\R peak detection and HRV\RAlps Rpeak detection');
SynchPath=(' C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Synchronizing data');
datacreationPath='C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Create data';
annotationfolder= ['C:\Users\310122653\Desktop\Fuckall of that shit\Annotation\participant' num2str(pat)];
annotationfiles=dir([annotationfolder '\*' num2str(epoch) '.csv ' ]);








dpath='E:\cECG_study\';
patient_folder_DAQ='20120809_test14';
DAQ_folder='20120809_130806';
Intelleviue_folder='2012-08-09 13-08_001';
% 
% patient_folder_DAQ='20120910_test17';
% DAQ_folder='20120910_122303';
% Intelleviue_folder='2012-09-10 12-23_001';

patient_folder_intellivue=[dpath patient_folder_DAQ '\Intelliview serial datalogging\'];

dataSetpath{4}= ([dpath patient_folder_DAQ '\']); % only patient
dataSetpath{1} = ([dpath patient_folder_DAQ '\' DAQ_folder]);
cd(dataSetpath{1})
dataSetpath{2} = ([patient_folder_intellivue Intelleviue_folder]); 






%% ***** READIN FILE *****
% Read in header
fileID=fopen([annotationfiles.folder '\' annotationfiles.name],'r');
tline = fgetl(fileID);
header(1,:) = regexp(tline, '\,', 'split');%  Split header

%Readin annotations (with first line header)
for ctr=2:9
    if ischar(tline)    
%           ctr = ctr + 1;
          tline = fgetl(fileID);         
          header(ctr,:) = regexp(tline, '\,', 'split'); 
    else
          break;     
    end
end

%  Parse and read rest of file
ctr = 10;
while(~feof(fileID))
    if ischar(tline)    
          ctr = ctr + 1;
          tline = fgetl(fileID);         
          annot(ctr-10,:) = regexp(tline, '\,', 'split'); 
    else
          break;     
    end
end

fclose(fileID);
clearvars ctr same tline 


% separate annotations per session
Sessionname=annot{2,1}; k=1;
Sessionnames{k}=annot{2,1}; %collecting all session names
for i=2:length(annot)
    same=strcmp(annot{i,1},Sessionname); % compare if the name is still the same
    if same==0 % if not...
        endsessions(k,1)=i; %write index in file...
        Sessionname=annot{i,1}; % and change Sessionname to the new filename
        Sessionnames{k}=Sessionname; %collecting all sessionames
        k=k+1;
    end
end
        
if exist('endsessions', 'var')==0 % if only one session exist, no endsession s created
    endsessions=length(annot);
end

%% ***** Creating sleep state array *****

% all annnotations for all sessions in one variable (wake)
wake= zeros(length(annot),1);
for i=1:length(annot)
    same=strcmp(annot{i,10},'		Wake');
  if same==1
      wake(i-1,1)=1;%-1 as first is header
  end
end

% all annnotations for all sessions in separate variables (Session_xxx)
for i=1:length(Sessionnames)
    if i==1 % for the first start with row 1
        Sessionnames{i}=regexprep(Sessionnames{i},{'\.','avi'},{''}); % delete dots from name      
        eval(['Session_' Sessionnames{i} '= wake(1:endsessions(i),1);' ])
    else
        eval(['Session_' Sessionnames{i} '= wake(endsessions(i-1):endsessions(i)+1,1);' ])        
    end
    
end

% all annnotations for all sessions in one variable (wake)
CareTaking= zeros(length(annot),1);
for i=1:length(annot)
    same=strcmp(annot{i,10},'		CareTaking');
  if same==1
      CareTaking(i-1,1)=1;%-1 as first is header
  end
end

% all annnotations for all sessions in separate variables (Session_xxx)
for i=1:length(Sessionnames)
    if i==1 % for the first start with row 1
        Sessionnames{i}=regexprep(Sessionnames{i},{'\.','avi'},{''}); % delete dots from name      
        eval(['Session_' Sessionnames{i} '= CareTaking(1:endsessions(i),1);' ])
    else
        eval(['Session_' Sessionnames{i} '= CareTaking(endsessions(i-1):endsessions(i)+1,1);' ])        
    end
    
end

%% FIND SHIFT IN ARRAY
%% load ECG

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
                disp([num2str(numFiles),' ECG bin files where found. File ' num2str(fileNumber) ' with a  lenght of ', num2str(ECG_bin_length) ])     
        end 
        %calculate FS from txt file
        cd(dataSetpath{1})
        fid = fopen('Ref ECG1_00000000.txt', 'r');
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
    
%     ECG_bin_tot_500=[zeros(900000,1); ECG_bin_tot_500];
    tECG=linspace(0,length(ECG_bin_tot_500)/500,length(ECG_bin_tot_500));

    clearvars ECG_bin_data_tot ECG_bin_cell

    
%% 

for i=1:length(Sessionnames)
    if i==1 % for the first start with row 1
        Sessionnames{i}=regexprep(Sessionnames{i},{'\.','avi'},{''}); % delete dots from name      
        eval(['Session_' Sessionnames{i} '= CareTaking(1:endsessions(i),1);' ])
        t{i,1}=linspace(0,eval(['length(Session_' Sessionnames{i} ')']),eval(['length(Session_' Sessionnames{i} ')']))';
    else
        eval(['Session_' Sessionnames{i} '= CareTaking(endsessions(i-1):endsessions(i)+1,1);' ])   
        t{i,1}=linspace(0,eval(['length(Session_' Sessionnames{i} ')']),eval(['length(Session_' Sessionnames{i} ')'])*500)';

    end
    
end

%% Plotting
adding=zeros(10,1);
eval(['actual_session=Session_' Sessionnames{i} ';']);
actual_session=actual_session*(0.6*max(ECG_bin_tot_500));
tECG=linspace(0,length(ECG_bin_tot_500)/500,length(ECG_bin_tot_500))';
Key_is_pressed=0;
%plotting

scrz=get(groot,'ScreenSize');
fh=figure('name','ASAB','Position',[20 scrz(4)/4 scrz(4)/0.6 scrz(4)/2]);
ah=axes('position',[0.05 0.05 0.9 0.9]); %left bottom width heigth
for j=1:length(ECG_bin_tot_500)
    plot(tECG,ECG_bin_tot_500); hold on
    xlim([0 5500])
    plot(t{1,1},actual_session);hold off
    xlim([0 5500])
    actual_session=[adding; actual_session(1:end-10,1)];
    title (['added sec:' num2str(j*10)]);
    pause(0.3)
end



