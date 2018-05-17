%%
%THIS  M FILE CREATE FEATURE MATRICES WHERE EVERYTHING IS COMBINED
% FOR MATRICES WHERE ECG, HR, BREATHING ETC ARE SEPARATE MATRICES USE Create_single_F_and_A_Matrix_Session_and_perPatient_forDNN
%%
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
tic
RRMethod='R'; %M or R if Michiels or Ralphs RR peak detection method was used 
dataset='MMC'; % ECG or cECG or MMC. In the fUture mayebe MMC and InnerSense
saving=1;
win=30; % window of annotations. 30 precicse, 300 smoothed
Datapack='Features';   %ECG_HRV_EDR, Features ,allRaw 


if strcmp('ECG',dataset) || strcmp('cECG',dataset)
    Pat=[4,5,6,7,9,10,11,12,13];
    datset='cECG_study\C_Processed_Data\';
    path='E:\';
    if strcmp('ECG',dataset)
        datapath=[path 'cECG_study\C_Processed_Data\HRV_features\'];
        datapathload=datset;
    elseif strcmp('cECG',dataset)
        datapath=[path 'cECG_study\C_Processed_Data\cHRV_features\'];    
        datapathload=datset;
    end
elseif strcmp('MMC',dataset) 
    PatientID=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]; % core. Show all patients in the folder
    Pat=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]; 


    datset='Users\310122653\Documents\PhD\Article_4_(MMC)\';
    path='C:\';
    datapathload=[ datset 'Processed data\'];
    datapath=[path datset 'Processed data\HRV_features\'];
end

        
if strcmp(RRMethod,'R')
    if     strcmp('ECG',dataset) || strcmp('MMC',dataset)
        
        savefolder= ([path datset 'DNN-Matrices\Matrices_' Datapack '\']);
        savefolderSession=([path datset 'DNN-Matrices\Matrices_' Datapack '\Sessions\']);    
      
    elseif strcmp('cECG',dataset)==1 && strcmp('Features',Datapack)
        savefolder= ([path datset '\DNN-Matrices\cMatrices_' Datapack '\']);
        savefolderSession=([path datset 'DNN-Matrices\cMatrices_' Datapack '\Sessions\']);    
    end

    TFeature_path=[datapath 'timedomain\'];
    FFeature_path=[datapath 'freqdomain\'];
    NLFeature_path=[datapath 'nonlinear\'];    

elseif strcmp(RRMethod,'M')
    if strcmp('ECG',dataset) || strcmp('MMC',dataset) 
        savefolder= ([path datset 'DNN-Matrices\MatricesM_' Datapack '\']);
        savefolderSession=([path datset 'DNN-Matrices\MatricesM_' Datapack '\Sessions\']);    
      
    elseif strcmp('cECG',dataset)==1 && strcmp('Features',Datapack)
        savefolder= ([path datset '\DNN-Matrices\cMatricesM_' Datapack '\']);
        savefolderSession=([path datset 'DNN-Matrices\cMatricesM_' Datapack '\Sessions\']);    
    end
    
    TFeature_path=[datapath 'timedomainM\'];
    FFeature_path=[datapath 'freqdomainM\'];
    NLFeature_path=[datapath 'nonlinearM\'];    
end

loadfolderAnnotation= [path datapathload 'Annotations\'];

if (exist(savefolder) )==0;  mkdir([savefolder]);end
if (exist(savefolderSession) )==0;  mkdir([savefolderSession]);end




if strcmp('allRaw',Datapack)==1
    Featurenames_time={...
        'ECG';...
        'HRV';...
        'EDR' ;...
        'Resp';...
        };
         Featurenames_frequency={};
         Featurenames_nonlinear={};
end

if strcmp('ECG_HRV_EDR',Datapack)==1
    Featurenames_time={...
        'ECG';...
        'HRV';...
        'EDR' ;...
        };
         Featurenames_frequency={};
         Featurenames_nonlinear={};
end         

if strcmp('Features',Datapack)==1
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
end
    

for N=1:length(Pat)
    disp(['Working on patient ' num2str(Pat(N))])
    disp('-------------------------------')
    Neonate=Pat(N);
%     N_I=find(PatientID==Neonate); % IF we do not start with 1 we have to choose the correct index
    Annottmp=[];FeatureMatrix_tmp=[]; % renew for each session. For MMC data we stack all sessions as we mostly only have 1
   
    FeatureMatrix={};tmp={};tmp2={};tmp3={};tmp4={};
    
    
    if strcmp('ECG',dataset) || strcmp('cECG',dataset)
        searchdiretory=['_Session_*_win_' num2str(win) '_*_'  num2str(Neonate) '.mat'];
        searchdiretoryT=['_win_' num2str(win) '_*_' num2str(Neonate) '.mat'];
        searchdiretoryF=['_win_' num2str(win) '_*_' num2str(Neonate) '.mat'];
        searchdirectoryNL=['_win_' num2str(win) '_*_' num2str(Neonate) '.mat'];
        searchdirectoryA1=['_win_' num2str(win) '_Intellivue_*_' num2str(Neonate) '.mat'];
    elseif strcmp('MMC',dataset) 
        searchdiretory=['_Session_*_win_' num2str(win) '_'  num2str(Neonate) '.mat'];
        searchdiretoryT=['_win_' num2str(win) '_' num2str(Neonate) '.mat'];
        searchdiretoryF=['_win_' num2str(win) '_' num2str(Neonate) '.mat'];
        searchdirectoryNL=['_win_' num2str(win) '_' num2str(Neonate) '.mat'];
        searchdirectoryA2=[loadfolderAnnotation 'Annotations_win_' num2str(win) '_' num2str(Neonate) '.mat'];    
    end

%--------------- Per Patient ---------------     
%#1 How many sessions?
%#2 Go through each session and load all Featrues and the annotation
%#3 Cut annotation and Feature to the same length
%#4 Merge Feature Sessions together or safe them as session

    Sessionnames=dir([TFeature_path Featurenames_time{1,1} searchdiretory]); 
    Sessionlength=length(cellfun('isempty',{Sessionnames.name}));

    for i=1:Sessionlength  
        disp(['Session ' num2str(i) '/' num2str(Sessionlength)])
        if strcmp('ECG',dataset) || strcmp('cECG',dataset)
            dateiname=dir([loadfolderAnnotation 'Annotations_Session_' num2str(i) searchdirectoryA1]);
        elseif strcmp('MMC',dataset) 
            dateiname=dir( searchdirectoryA2);
        end
        load([loadfolderAnnotation dateiname.name]);


        
    % all from one patient TIME DOMAIN
        for j=1:length(Featurenames_time) 
            dateiname=dir([TFeature_path Featurenames_time{j,1} '_Session_' num2str(i) searchdiretoryT]);
            load([TFeature_path dateiname.name])
            if iscolumn(Feature)
                Feature=Feature';% Need to be a row vector
            end             
            if strcmp('MMC',dataset) && Neonate==17
                 Feature=Feature(1,2:end); % as in pat 17 of the MMC data the annotations start exactly 30s after the ECG, we first cut the first cell of the feature to have the same start
            end
            if length(Feature)>length(Annotations)
                Feature=Feature(1:length(Annotations)); % Cut the Feature to the length of the annotation. We asume that the annotations always start at the beginning but end earlier
            elseif length(Annotations)>length(Feature)
                Annotations=Annotations(1:length(Feature));
            end

            for f=1:length(Feature)
                if iscell(Feature)==1 && isrow(Feature{1,f})
                    Feature{1,f}=Feature{1,f}';
                end
            end
            if iscell(Feature)==1 % Some features are saved in cell some as double
%                 tmp= {tmp; cell2mat(Feature)}; % adding each session of one feature to one single line of Feature 
                tmp= [tmp; Feature]; % adding each session of one feature to one single line of Feature 
                
            else
                tmp= [tmp; (num2cell(Feature))]; % adding each session of one feature to one single line of Feature 
            end             
        end
        tmp2=[tmp2,tmp]; % Adding the long single line of Feature into the Feature Matrix, where each row is one Feature        
        tmp={};

        % all from on patient FREQUENCY DOMAIN
        if strcmp('Features',Datapack) % the other datapacks do not have nonlinear and frequency and therefor problems arise
            for j=1:length(Featurenames_frequency) 
                dateiname=dir([FFeature_path Featurenames_frequency{j,1} '_Session_' num2str(i) searchdiretoryF]);
                load([FFeature_path dateiname.name]);
                if iscolumn(Feature) % Need to be a row vector
                    Feature=Feature';
                end                  
                if strcmp('MMC',dataset) && Neonate==17
                    Feature=Feature(1,2:end); % as in pat 17 of the MMC data the annotations start exactly 30s after the ECG, we first cut the first cell of the feature to have the same start
                end                
                if length(Feature)>length(Annotations)
                    Feature=Feature(1:length(Annotations)); % Cut the Feature to the length of the annotation. We asume that the annotations always start at the beginning but end earlier
                elseif length(Annotations)>length(Feature)
                    Annotations=Annotations(1:length(Feature));
                end
       
                for f=1:length(Feature)
                    if iscell(Feature)==1 && isrow(Feature{1,f})
                        Feature{1,f}=Feature{1,f}';
                    end
                end
                if iscell(Feature)==1 % Some features are saved in cell some as double
    %               tmp= {tmp; cell2mat(Feature)}; % adding each session of one feature to one single line of Feature 
                    tmp= [tmp; Feature]; % adding each session of one feature to one single line of Feature 
                else
                    tmp= [tmp; num2cell(Feature)]; % adding each session of one feature to one single line of Feature 
                end             
            end
            tmp3=[tmp3 ,tmp]; % Adding the long single line of Feature into the Feature Matrix, where each row is one Feature
            tmp={};
    
    
% all from on patient NONELINEAR    
            for j=1:length(Featurenames_nonlinear) 
                dateiname=dir([NLFeature_path Featurenames_nonlinear{j,1} '_Session_' num2str(i)  searchdirectoryNL]);
                load([NLFeature_path dateiname.name])
                if iscolumn(Feature) % Need to be a row vector
                    Feature=Feature';
                end
                if strcmp('MMC',dataset) && Neonate==17
                    Feature=Feature(1,2:end); % as in pat 17 of the MMC data the annotations start exactly 30s after the ECG, we first cut the first cell of the feature to have the same start
                end                
                if length(Feature)>length(Annotations)
                    Feature=Feature(1:length(Annotations)); % Cut the Feature to the length of the annotation. We asume that the annotations always start at the beginning but end earlier
                elseif length(Annotations)>length(Feature)
                    Annotations=Annotations(1:length(Feature));
                end

                for f=1:length(Feature) 
                    if iscell(Feature)==1 && isrow(Feature{1,f})
                        Feature{1,f}=Feature{1,f}';
                    end
                end            
                if iscell(Feature)==1 % Some features are saved in cell some as double
    %               tmp= {tmp; cell2mat(Feature)}; % adding each session of one feature to one single line of Feature 
                    tmp= [tmp; Feature]; % adding each session of one feature to one single line of Feature             
                else
                    tmp= [tmp; num2cell(Feature)]; % adding each session of one feature to one single line of Feature 
                end             
            end
            tmp4=[tmp4 ,tmp]; % Adding the long single line of Feature into the Feature Matrix, where each row is one Feature
            tmp={};
        end
        FeatureMatrix=[tmp2; tmp3 ;tmp4]; % Adding the long single line of Feature into the Feature Matrix, where each row is one Feature
        tmp2={};tmp3={};tmp4={}; % Resetting 

        % remove total epoch if all values are nan, also in annotations to
        % not lose synch
        idx=cellfun(@(FeatureMatrix) all(isnan(FeatureMatrix)),FeatureMatrix);
        c=any(idx);% find index where one cell element(any) is all nan (row before)
        FeatureMatrix(:,c)=[];
        Annotations(:,c)=[];
        clearvars c
        % remove all nans from each cell, as there are mostly less nans than 30s, it just
        % reduces the first epoch. Deleting nans as the standart scaler of Python cannot handle nans
        for fm = 1:numel(FeatureMatrix)
            FeatureMatrix{fm} = FeatureMatrix{fm}(~isnan(FeatureMatrix{fm}),:) ;            
        end
        FMtmp=cell2mat(FeatureMatrix);
        
        %STANDARD SCALE WITH Z-SCORE
        for zs=1:size(FMtmp,1)
            FMtmp(zs,:) = zscore(FMtmp(zs,:));
        end

        %interpolate to gain same length per Feature
%         FMtmp = cell(size(FeatureMatrix));
        
%         for k = 1:numel(FeatureMatrix)
%            if numel(FeatureMatrix{k})~=1 % interpolation needs at least two values
%                N = numel(FeatureMatrix{k});
%                xo = 1:N;
%                xi = linspace(1,N,win*500);
%                FMtmp{k} = interp1(xo,FeatureMatrix{k},xi);
%            else % if only one value (e.g. detected R peak ) due to noise, just create cell with this one value in the specific length 
%                FMtmp{k}=kron(FeatureMatrix{k}, ones(1,win*500));
%            end     
%         end


            %standart scale with (x-mean)/std 
%         for F=1:size(FMtmp,1)
%             flattenedM=cell2mat(FMtmp(F,:));
%             MeanStd(F,1)=mean(flattenedM);
%             MeanStd(F,2)=std(flattenedM);            
%         end
        % Cell to 3D array ---lets do that in Keras
        % IN THE FOLLOWING LINE YOU CAN DETERMINE THE ORDER OF THE DATA.
        % KERAS NEEDS [SAMPLES, TIMESTEPS, FEATURES]
%         FMtmp=permute(reshape(FMtmp.',numel(FMtmp(1,:)),size(FMtmp,1),[]),[3,1,2]);
%         FMtmp=permute(reshape(cell2mat(FMtmp).',numel(FMtmp{1}),size(FMtmp,2),[]),[1,3,2]);

        FeatureMatrix=FMtmp;

        if saving
            Saving_Session_F(FeatureMatrix,savefolderSession, Neonate,i)
            Saving_Session_A(Annotations,savefolderSession, Neonate,i)
        end
        
%         if strcmp('MMC',dataset) 
%             longestAnnotation=623;
%         end
%         
% % here zeros padd to create 1 3D tensor % NOOO DO ZERO PADDING IN PYTHON
% AS WE STILL HAVE TO BE VARIABLE HOW TO STACK AND SEPARATE THE DATA
%         if F<longestAnnotation % 623 is the longest data for the MMC data. F is the length of the Feature after removing nans
%             FeatureMatrix = padarray(FeatureMatrix,[0 623-F 0],'post') ;
%             Annotations = padarray(cell2mat(Annotations),[0 623-F],'post') ;
%         end
% 
            FeatureMatrix_tmp=[FeatureMatrix_tmp FeatureMatrix];
            Annottmp=[Annottmp Annotations]; % 
% 


    end % Session  
    
    Annotations=Annottmp;
    FeatureMatrix=FeatureMatrix_tmp;

    if saving
        Saving_F(FeatureMatrix,savefolder, Neonate,win)
        Saving_A(Annotations,savefolder, Neonate,win)
    end
    
    FeatureMatrix=[];
    Annotations=[];
end
disp('* * * * * * * * * * ')
toc



%% Nested saving
    function Saving_F(FeatureMatrix,savefolder, Neonate,win)
        if exist('FeatureMatrix','var')==1
            name=inputname(1); % variable name of function input
            save([savefolder name '_' num2str(Neonate)],'FeatureMatrix')
        else
            disp(['saving of ' name ' not possible'])
        end       
    end
    
    function Saving_Session_F(FeatureMatrix,savefolder, Neonate,Session,win)
        if exist('FeatureMatrix','var')==1
            name=inputname(1); % variable name of function input
            save([savefolder name '_' num2str(Neonate) '_Session_' num2str(Session) ],'FeatureMatrix')
        else
            disp(['saving of ' name ' not possible'])
        end       
    end
    
    function Saving_A(Annotations,savefolder, Neonate,win)
        if exist('Annotations','var')==1
            name=inputname(1); % variable name of function input
            save([savefolder name '_' num2str(Neonate)],'Annotations')
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