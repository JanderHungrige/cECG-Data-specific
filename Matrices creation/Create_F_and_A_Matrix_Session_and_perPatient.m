%18.4.2018 LAst change , stuff for Xi
%5.7.2018 last change, added z score per feature and per session 
%This m file generates FEature and annotation Matrices. Per patient and per
%session
% Each Row is one Feature. Where each session of the same Feature


%      '	1=	ActiveSleep';...
%         '	2=	QuietSleep';...
%         '	3=	Wake';...
%         '	4=	CareTaking';...
%         '	5=	UnknownBedState'...
%           6=  Transition


%#1 How many sessions?
%#2 Go through each session and load all Featrues and the annotation
%#3 Cut annotation and Feature to the same length
%#4 Save each combination as Session in multile matrices
%#5 Merge Feature Sessions together Safe them as one matrix
clear
clc


%1 UMBENENEN hrv FEATURES IN INTELLIVUE FEATURES
% run annotaions for 4 9 10 12 separate
% create matrice
% do same for DAQ (not 4)



RRMethod='R'; %M or R if Michiels or Ralphs RR peak detection method was used 
dataset='ECG'; % ECG or cECG. In the fUture mayebe MMC and InnerSense
saving=1;
win=300; % window of annotations. 30 precicse, 300 smoothed
win=30;
for_DNN=1 % If 1 then Annotations ara changed into cells to be th same as the once from MMC


Pat=[4,5,6,7,9,10,11,12,13];
% Pat=6;
Pat=[4,5,6,7,11,13];
% path='E:\';
Pat=[4,5,6,7,11,13];
Pat=[9,10,11,12,13]
Pat=[4,5,6,7,9,10,11,12,13];
Device='Intellivue' ;%DAQ Intellivue or something else
if strcmp(Device,'DAQ')
    datarecording='HRV_Features_DAQ';
    deviceMatric = 'Matrix_DAQ';
    disp('using DAQ')
elseif strcmp(Device,'Intellivue')
    datarecording='HRV_Features_Intellivue';
    deviceMatric = 'Matrix_Intellivue';
    disp('using Intellivue')
else
    datarecording='HRV_Features';
    deviceMatric = 'Matrix'; 
    disp('using overall')
end
path=('C:\Users\310122653\Documents\PhD\');
if strcmp('ECG',dataset)==1
%     datapath=[path 'cECG_study\C_Processed_Data\HRV_features\'];
    datapath=[path 'Article_3_(cECG)\Processed Data\'];
elseif strcmp('cECG',dataset)==1
%     datapath=[path 'cECG_study\C_Processed_Data\cHRV_features\']; 
    datapath=[path 'Article_3_(cECG)\Processed Data\'];
end
if strcmp(RRMethod,'R')
    if strcmp('ECG',dataset)==1
%         savefolder= ([path 'cECG_study\C_Processed_Data\Matrices\']);
%         savefolderSession=([path 'cECG_study\C_Processed_Data\Matrices\Sessions\']);    
        savefolder= ([datapath deviceMatric '\']); %inellivue
        savefolderSession=([datapath deviceMatric '\Sessions\']);    %inellivue
        mkdir(savefolder);
        mkdir(savefolderSession);
    elseif strcmp('cECG',dataset)==1
%         savefolder= ([path 'cECG_study\C_Processed_Data\cMatrices\']);
%         savefolderSession=([path 'cECG_study\C_Processed_Data\cMatrices\Sessions\']);    
        savefolder= ([datapath 'cMatrices\']);
        savefolderSession=([datapath 'cMatrices\Sessions\']);   
    end
%     TFeature_path=[datapath 'timedomain\'];
%     FFeature_path=[datapath 'freqdomain\'];
%     NLFeature_path=[datapath 'nonlinear\'];    
    TFeature_path=[datapath datarecording '\timedomain\'];
    FFeature_path=[datapath datarecording '\freqdomain\'];
    NLFeature_path=[datapath datarecording '\nonlinear\'];  
% elseif strcmp(RRMethod,'M')
%      if strcmp('ECG',dataset)==1
%         savefolder= ([path 'cECG_study\C_Processed_Data\MatricesM\']);
%         savefolderSession=([path 'cECG_study\C_Processed_Data\MatricesM\Sessions\']);    
% 
%     elseif strcmp('cECG',dataset)==1
%         savefolder= ([path 'cECG_study\C_Processed_Data\cMatricesM\']);
%         savefolderSession=([path 'cECG_study\C_Processed_Data\cMatricesM\Sessions\']);    
%      end 
%     TFeature_path=[datapath 'timedomainM\'];
%     FFeature_path=[datapath 'freqdomainM\'];
%     NLFeature_path=[datapath 'nonlinearM\'];
    
end


% loadfolderAnnotation= [path 'cECG_study\C_Processed_Data\Annotations\'];
loadfolderAnnotation= [datapath 'Annotations\'];


Featurenames_time={...
    'BpE';...
    'lineLength';...
    'meanlineLength';...
    'NN10'; 'NN20';'NN30';'NN50';...
    'pNN10'; 'pNN20';'pNN30';'pNN50';...
    'RMSSD';...
    'SDaLL';...
    'SDANN';...
    'SDLL';...
    'SDNN';...
    'pDEC';...
    'SDDEC';...    
    };

Featurenames_frequency={...
    'HF';...
    'HFnorm';...
    'LF';...
    'LFnorm';...
    'ratioLFHF';...
    'sHF';...
    'sHFnorm';...
    'totpow';...
    'uHF';...
    'uHFnorm';...
    'VLF';...
    };

Featurenames_nonlinear={...
    'SampEn';...
    'QSE';...
    'SEAUC';...    
    'LZNN';...
    'LZECG';... 
    };

for N=1:length(Pat)
    disp(['Working on patient ' num2str(Pat(N))])
    disp('-------------------------------')
    Neonate=Pat(N);
    Annottmp=[];FeatureMatrix_tmp=[];
    FeatureMatrix=[];tmp=[];tmp2=[];tmp3=[];tmp4=[];
%--------------- Per Patient ---------------     
%#1 How many sessions?
%#2 Go through each session and load all Featrues and the annotation
%#3 Cut annotation and Feature to the same length
%#4 Merge Feature Sessions together or safe them as session

Sessionlength=dir([TFeature_path Featurenames_time{1,1} '_Session_*_' num2str(Neonate) '.mat']); 
Sessionlength=length(cellfun('isempty',{Sessionlength.name}));

    for i=1:Sessionlength  
%         dateiname=dir([loadfolderAnnotation 'Annotations_Session_' num2str(i) '_win_' num2str(win) '_Intellivue_*_' num2str(Neonate) '.mat']);
        dateiname=dir([loadfolderAnnotation 'Annotations_Session_' num2str(i) '_pat_' num2str(Neonate) '.mat']);

        load([loadfolderAnnotation dateiname.name]);
    % all from one patient TIME DOMAIN
        for j=1:length(Featurenames_time) 
%             dateiname=dir([TFeature_path Featurenames_time{j,1} '_Session_' num2str(i) '_*_' num2str(Neonate) '.mat']);
            dateiname=dir([TFeature_path Featurenames_time{j,1} '_Session_' num2str(i) '_*_' num2str(Neonate) '.mat']);

            load([TFeature_path dateiname.name])
            if length(Feature)>length(Annotations)
                Feature=Feature(1:length(Annotations)); % Cut the Feature to the length of the annotation. We asume that the annotations always start at the beginning but end earlier
            elseif length(Annotations)>length(Feature)
                Annotations=Annotations(1:length(Feature));
            end
            if iscolumn(Feature)
                Feature=Feature';% Need to be a row vector
            end              
            if iscell(Feature)==1 % Some features are saved in cell some as double
                tmp= [tmp; cell2mat(Feature)]; % adding each session of one feature to one single line of Feature 
            else
                tmp= [tmp; Feature]; % adding each session of one feature to one single line of Feature 
            end 
        end
        tmp=zscore(tmp,1,2);% per sample(1), per row(2) Standardizing the data to mean 0 and standard deviation                
        tmp2=[tmp2 ,tmp]; % Adding the long single line of Feature into the Feature Matrix, where each row is one Feature        
        tmp=[];

    % all from on patient FREQUENCY DOMAIN
        for j=1:length(Featurenames_frequency) 
%             dateiname=dir([FFeature_path Featurenames_frequency{j,1} '_Session_' num2str(i) '_*_' num2str(Neonate) '.mat']);
            dateiname=dir([FFeature_path Featurenames_frequency{j,1} '_Session_' num2str(i) '_*_' num2str(Neonate) '.mat']);
            
            load([FFeature_path dateiname.name]);
            if length(Feature)>length(Annotations)
                Feature=Feature(1:length(Annotations)); % Cut the Feature to the length of the annotation. We asume that the annotations always start at the beginning but end earlier
            elseif length(Annotations)>length(Feature)
                Annotations=Annotations(1:length(Feature));
            end
            if iscolumn(Feature) % Need to be a row vector
                Feature=Feature';
            end              
            if iscell(Feature)==1 % Some features are saved in cell some as double
            tmp= [tmp; cell2mat(Feature)]; % adding each session of one feature to one single line of Feature 
            else
            tmp= [tmp; Feature]; % adding each session of one feature to one single line of Feature 
            end  
        end
        tmp=zscore(tmp,1,2);% per sample(1), per row(2) Standardizing the data to mean 0 and standard deviation                
        tmp3=[tmp3 ,tmp]; % Adding the long single line of Feature into the Feature Matrix, where each row is one Feature
        tmp=[];
    
    
% all from on patient NONELINEAR    
        for j=1:length(Featurenames_nonlinear) 
%             dateiname=dir([NLFeature_path Featurenames_nonlinear{j,1} '_Session_' num2str(i) '_*_' num2str(Neonate) '.mat']);
            dateiname=dir([NLFeature_path Featurenames_nonlinear{j,1} '_Session_' num2str(i) '_*_' num2str(Neonate) '.mat']);                        
            load([NLFeature_path dateiname.name])
            if length(Feature)>length(Annotations)
                Feature=Feature(1:length(Annotations)); % Cut the Feature to the length of the annotation. We asume that the annotations always start at the beginning but end earlier
            elseif length(Annotations)>length(Feature)
                Annotations=Annotations(1:length(Feature));
            end
            if iscolumn(Feature) % Need to be a row vector
                Feature=Feature';
            end                
            if iscell(Feature)==1 % Some features are saved in cell some as double
            tmp= [tmp; cell2mat(Feature)]; % adding each session of one feature to one single line of Feature 
            else
            tmp= [tmp; Feature]; % adding each session of one feature to one single line of Feature 
            end 
        end
        tmp=zscore(tmp,1,2);% per sample(1), per row(2) Standardizing the data to mean 0 and standard deviation                
        tmp4=[tmp4 ,tmp]; % Adding the long single line of Feature into the Feature Matrix, where each row is one Feature
        tmp=[];
        FeatureMatrix=[tmp2; tmp3 ;tmp4]; % Adding the long single line of Feature into the Feature Matrix, where each row is one Feature
        tmp2=[];tmp3=[];tmp4=[]; % Resetting 
        
        % Deleting nans as the standart scaler of Python cannot handle nans
%         Simultaniousely in Features and Annotations to not loose synch
        [r,c]=find(any(isnan(FeatureMatrix),1));
        Annotations=cell2mat(Annotations);
        FeatureMatrix(:,c)=[];
        Annotations(:,c)=[];
        clearvars r c
        

        if saving
            Saving_Session_F(FeatureMatrix,savefolderSession, Neonate,i)
            Saving_Session_A(Annotations,savefolderSession, Neonate,i)
        end
        
        FeatureMatrix_tmp=[FeatureMatrix_tmp FeatureMatrix];
        Annottmp=[Annottmp Annotations]; % 

    end % Session       
    Annotations=Annottmp; 
    FeatureMatrix=FeatureMatrix_tmp;
    % The FEature Matrix is composed of standardised sessions. Therfore
    % mean is laready 0 and std 1 . Another z-score does not alter the
    % values anymore. We leave it in that it is not raised again.
%     FeatureMatrix=zscore(FeatureMatrix,1,2);% per sample(1), per row(2) Standardizing the data to mean 0 and standard deviation                
if saving
    Saving_F(FeatureMatrix,savefolder, Neonate)
    Saving_A(Annotations,savefolder, Neonate)
    
end
disp(length(FeatureMatrix))

FeatureMatrix=[];
Annotations=[];
end


%% Nested saving
    function Saving_F(FeatureMatrix,savefolder, Neonate)
        if exist('FeatureMatrix','var')==1
            name=inputname(1); % variable name of function input
            save([savefolder name '_' num2str(Neonate) '_win_30'],'FeatureMatrix')
        else
            disp(['saving of ' name ' not possible'])
        end       
    end
    
    function Saving_Session_F(FeatureMatrix,savefolder, Neonate,Session)
        if exist('FeatureMatrix','var')==1
            name=inputname(1); % variable name of function input
            save([savefolder name '_' num2str(Neonate) '_Session_' num2str(Session) ],'FeatureMatrix')
        else
            disp(['saving of ' name ' not possible'])
        end       
    end
    
    function Saving_A(Annotations,savefolder, Neonate)
        if exist('Annotations','var')==1
            name=inputname(1); % variable name of function input
            save([savefolder name '_' num2str(Neonate) '_win_30'],'Annotations')
        else
            disp(['saving of ' name ' not possible'])
        end       
    end    
    
    function Saving_Session_A(Annotations,savefolder, Neonate,Session)
        if exist('Annotations','var')==1
            name=inputname(1); % variable name of function input
            save([savefolder name '_' num2str(Neonate) '_Session_' num2str(Session) ],'Annotations')
        else
            disp(['saving of ' name ' not possible'])
        end       
    end
    



% end