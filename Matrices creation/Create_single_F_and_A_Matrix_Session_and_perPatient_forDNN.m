%%
%THIS  M FILE CREATE FEATURE MATRICES WHERE ECG, HR, BREATHING ETC ARE SEPARATE MATRICES
%% HEADER
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
%% DECLARATIONS
RRMethod='R'; %M or R if Michiels or Ralphs RR peak detection method was used 
dataset='ECG'; % ECG or cECG. In the fUture mayebe MMC and InnerSense
saving=1;
Datapack=['ECG';'HRV';'EDR';'RSP';'FET'] ;  
Pat=[4,5,6,7,9,10,11,12,13];

path='E:\';
if strcmp('ECG',dataset)==1
    datapath=[path 'cECG_study\C_Processed_Data\HRV_features\'];
elseif strcmp('cECG',dataset)==1
    datapath=[path 'cECG_study\C_Processed_Data\cHRV_features\'];    
end
TFeature_path=[datapath 'timedomain\'];
FFeature_path=[datapath 'freqdomain\'];
NLFeature_path=[datapath 'nonlinear\'];    
loadfolderAnnotation= [path 'cECG_study\C_Processed_Data\Annotations\'];


    Featurenames_time={...
        'ECG';...
        'HRV';...
        'EDR' ;...
        'Resp';...
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

for DT=1:numel(Datapack)/3
    if strcmp(Datapack(DT,:),'ECG')==1 
        win=30; % window of annotations. 30 precicse, 300 smoothed
        Featurenames_time_use{1}=Featurenames_time{1};
        Featurenames_frequency_use={};
        Featurenames_nonlinear_use={};
        if strcmp('ECG',dataset)==1
            savefolder=([path 'cECG_study\C_Processed_Data\DNN-Matrices\Matrices\' Datapack(DT,:) '\' ]);
            savefolderSession=([path 'cECG_study\C_Processed_Data\DNN-Matrices\Matrices\' Datapack(DT,:) '\Sessions\']);
        elseif strcmp('cECG',dataset)==1
            savefolder=([path 'cECG_study\C_Processed_Data\DNN-Matrices\cMatrices\' Datapack(DT,:) '\']);
            savefolderSession=([path 'cECG_study\C_Processed_Data\DNN-Matrices\cMatrices\' Datapack(DT,:) '\Sessions\']);  
        end
    elseif strcmp(Datapack(DT,:),'HRV')==1
        win=30; % window of annotations. 30 precicse, 300 smoothed
        Featurenames_time_use{1}=Featurenames_time{2};
        Featurenames_frequency_use={};
        Featurenames_nonlinear_use={}; 
        if strcmp('ECG',dataset)==1
            savefolder=([path 'cECG_study\C_Processed_Data\DNN-Matrices\Matrices\' Datapack(DT,:) '\' ]);
            savefolderSession=([path 'cECG_study\C_Processed_Data\DNN-Matrices\Matrices\' Datapack(DT,:) '\Sessions\']);
        elseif strcmp('cECG',dataset)==1
            savefolder=([path 'cECG_study\C_Processed_Data\DNN-Matrices\cMatrices\' Datapack(DT,:) '\']);
            savefolderSession=([path 'cECG_study\C_Processed_Data\DNN-Matrices\cMatrices\' Datapack(DT,:) '\Sessions\']);    
        end         
    elseif strcmp(Datapack(DT,:),'EDR')==1
        win=30; % window of annotations. 30 precicse, 300 smoothed
        Featurenames_time_use{1}=Featurenames_time{3};
        Featurenames_frequency_use={};
        Featurenames_nonlinear_use={}; 
        if strcmp('ECG',dataset)==1
            savefolder=([path 'cECG_study\C_Processed_Data\DNN-Matrices\Matrices\' Datapack(DT,:) '\' ]);
            savefolderSession=([path 'cECG_study\C_Processed_Data\DNN-Matrices\Matrices\' Datapack(DT,:) '\Sessions\']);
        elseif strcmp('cECG',dataset)==1
            savefolder=([path 'cECG_study\C_Processed_Data\DNN-Matrices\cMatrices\' Datapack(DT,:) '\']);
            savefolderSession=([path 'cECG_study\C_Processed_Data\DNN-Matrices\cMatrices\' Datapack(DT,:) '\Sessions\']);   
        end         
    elseif strcmp(Datapack(DT,:),'RSP')==1 
        win=30; % window of annotations. 30 precicse, 300 smoothed
        Featurenames_time_use{1}=Featurenames_time{4};
        Featurenames_frequency_use={};
        Featurenames_nonlinear_use={}; 
        if strcmp('ECG',dataset)==1
            savefolder=([path 'cECG_study\C_Processed_Data\DNN-Matrices\Matrices\' Datapack(DT,:) '\' ]);
            savefolderSession=([path 'cECG_study\C_Processed_Data\DNN-Matrices\Matrices\' Datapack(DT,:) '\Sessions\']);
        elseif strcmp('cECG',dataset)==1
            savefolder=([path 'cECG_study\C_Processed_Data\DNN-Matrices\cMatrices\' Datapack(DT,:) '\']);
            savefolderSession=([path 'cECG_study\C_Processed_Data\DNN-Matrices\cMatrices\' Datapack(DT,:) '\Sessions\']);  
        end  
    elseif strcmp(Datapack(DT,:),'PWR')==1 
        win=30; % window of annotations. 30 precicse, 300 smoothed
        Featurenames_time_use{1}=Featurenames_time{4};
        Featurenames_frequency_use={};
        Featurenames_nonlinear_use={}; 
        if strcmp('ECG',dataset)==1
            savefolder=([path 'cECG_study\C_Processed_Data\DNN-Matrices\Matrices\' Datapack(DT,:) '\' ]);
            savefolderSession=([path 'cECG_study\C_Processed_Data\DNN-Matrices\Matrices\' Datapack(DT,:) '\Sessions\']);
        elseif strcmp('cECG',dataset)==1
            savefolder=([path 'cECG_study\C_Processed_Data\DNN-Matrices\cMatrices\' Datapack(DT,:) '\']);
            savefolderSession=([path 'cECG_study\C_Processed_Data\DNN-Matrices\cMatrices\' Datapack(DT,:) '\Sessions\']);  
        end         
    elseif strcmp(Datapack(DT,:),'FET')==1 
        win=300;
        Featurenames_time_use=Featurenames_time(5:end);
        Featurenames_frequency_use=Featurenames_frequency;
        Featurenames_nonlinear_use=Featurenames_nonlinear; 
        if strcmp('ECG',dataset)==1
            savefolder=([path 'cECG_study\C_Processed_Data\DNN-Matrices\Matrices\Datapack(DT,:)']);
            savefolder=([path 'cECG_study\C_Processed_Data\DNN-Matrices\Matrices\' Datapack(DT,:) '\' ]);
            savefolderSession=([path 'cECG_study\C_Processed_Data\DNN-Matrices\Matrices\' Datapack(DT,:) '\Sessions\']);
        elseif strcmp('cECG',dataset)==1
            savefolder=([path 'cECG_study\C_Processed_Data\DNN-Matrices\cMatrices\' Datapack(DT,:) '\']);
            savefolderSession=([path 'cECG_study\C_Processed_Data\DNN-Matrices\cMatrices\' Datapack(DT,:) '\Sessions\']);    
        end         
    end
%% PROCESSING  
    for N=1:length(Pat)
        disp('-------------------------------')        
        disp(['Working on patient ' num2str(Pat(N))])        
        disp(Datapack(DT,:))        
        Neonate=Pat(N);
        Annottmp=[];FeatureMatrix_tmp=[];
        FeatureMatrix={};tmp={};tmp2={};tmp3={};tmp4={};
    %--------------- Per Patient ---------------     
    %#1 How many sessions?
    %#2 Go through each session and load all Featrues and the annotation
    %#3 Cut annotation and Feature to the same length
    %#4 Merge Feature Sessions together or safe them as session

    Sessionlength=dir([TFeature_path Featurenames_time_use{1,1} '_Session_*_win_' num2str(win) '_*_'  num2str(Neonate) '.mat']); 
    Sessionlength=length(cellfun('isempty',{Sessionlength.name}));

        for i=1:Sessionlength  
            disp(['Session ' num2str(i) '/' num2str(Sessionlength)])

            dateiname=dir([loadfolderAnnotation 'Annotations_Session_' num2str(i) '_win_' num2str(win) '_Intellivue_*_' num2str(Neonate) '.mat']);
            load([loadfolderAnnotation dateiname.name]);

        % all from one patient TIME DOMAIN
            for j=1:length(Featurenames_time_use) 
                dateiname=dir([TFeature_path Featurenames_time_use{j,1} '_Session_' num2str(i) '_win_' num2str(win) '_*_' num2str(Neonate) '.mat']);
                load([TFeature_path dateiname.name])
                if length(Feature)>length(Annotations)
                    Feature=Feature(1:length(Annotations)); % Cut the Feature to the length of the annotation. We asume that the annotations always start at the beginning but end earlier
                elseif length(Annotations)>length(Feature)
                    Annotations=Annotations(1:length(Feature));
                end
                if iscolumn(Feature)
                    Feature=Feature';% Need to be a row vector
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
        if strcmp('FET',Datapack(DT,:))==1  
            for j=1:length(Featurenames_frequency_use) 
                dateiname=dir([FFeature_path Featurenames_frequency_use{j,1} '_Session_' num2str(i) '_win_' num2str(win) '_*_' num2str(Neonate) '.mat']);
                load([FFeature_path dateiname.name]);
                if length(Feature)>length(Annotations)
                    Feature=Feature(1:length(Annotations)); % Cut the Feature to the length of the annotation. We asume that the annotations always start at the beginning but end earlier
                elseif length(Annotations)>length(Feature)
                    Annotations=Annotations(1:length(Feature));
                end
                if iscolumn(Feature) % Need to be a row vector
                    Feature=Feature';
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
                dateiname=dir([NLFeature_path Featurenames_nonlinear{j,1} '_Session_' num2str(i) '_win_' num2str(win) '_*_' num2str(Neonate) '.mat']);
                load([NLFeature_path dateiname.name])
                if length(Feature)>length(Annotations)
                    Feature=Feature(1:length(Annotations)); % Cut the Feature to the length of the annotation. We asume that the annotations always start at the beginning but end earlier
                elseif length(Annotations)>length(Feature)
                    Annotations=Annotations(1:length(Feature));
                end
                if iscolumn(Feature) % Need to be a row vector
                    Feature=Feature';
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
        end %strcmp FET            
        
            FeatureMatrix=[tmp2; tmp3 ;tmp4]; % Adding the long single line of Feature into the Feature Matrix, where each row is one Feature
            tmp2={};tmp3={};tmp4={}; % Resetting 

            % remove total epoch if all values are nan, also in annotations to
            % not lose synch
%             c=any(idx);% find index where one cell element(any) is all nan (row before)
            if size(FeatureMatrix,1)==1 %if FeatureMatrix = 1xX Matrix
                idx=cellfun(@(FeatureMatrix) all(isnan(FeatureMatrix)),FeatureMatrix);                
            else % if FeatureMatrix = XxY Matrix
                idx=find(any(isnan(cell2mat(FeatureMatrix)),1));                              
            end
            FeatureMatrix(:,idx)=[];
            Annotations(:,idx)=[];              
            clearvars idx
            % remove all nans from each cell, as there are mostly less nans than 30s, it just
            % reduces the first epoch. Deleting nans as the standart scaler of Python cannot handle nans
            for fm = 1:numel(FeatureMatrix)
                FeatureMatrix{fm} = FeatureMatrix{fm}(~isnan(FeatureMatrix{fm}),:) ;            
            end

            %interpolate to gain same length per Feature
            if strcmp(Datapack(DT,:),'FET')~=1
                FMtmp = cell(size(FeatureMatrix));
                for k = 1:numel(FeatureMatrix)
                   if numel(FeatureMatrix{k})~=1 % interpolation needs at least two values
                       N = numel(FeatureMatrix{k});
                       xo = 1:N;
                       xi = linspace(1,N,win*500);
                       FMtmp{k} = interp1(xo,FeatureMatrix{k},xi)' ;
                   else % if only one value (e.g. detected R peak ) due to noise, just create cell with this one value in the specific length 
                       FMtmp{k}=kron(FeatureMatrix{k}, ones(1,win*500))' ;
                   end     
                end
            elseif strcmp(Datapack(DT,:),'FET')==1 % do not interpolate the single features
                FMtmp=FeatureMatrix;  
            end
            
            
    % normalize to no mean std 1 with (x-mean)/std 
            if size(FMtmp,1)==1 % only one feature e.g. ECG                
                flattenedM=cell2mat(FMtmp);
                MeanStd(1,1)=mean2(flattenedM);
                MeanStd(1,2)=std2(flattenedM);   
            else
                for F=1:size(FMtmp,1) % multible features 
                    MeanStd(F,1)=nanmean(cell2mat(FMtmp(F,:)));
                    MeanStd(F,2)=nanstd(cell2mat(FMtmp(F,:)));
                end                
            end
            % Cell to 3D array
            % With permut you can transpose a 3D Matrix.
            % KERAS NEEDS [SAMPLES, TIMESTEPS, FEATURES]
            if strcmp(Datapack(DT,:),'FET')~=1 % all else than FET
                FMtmp=permute(reshape(cell2mat(FMtmp).',numel(FMtmp{1}),size(FMtmp,2),[]),[2,1,3]);
                for Samp=1:size(FMtmp,1) % used with [2 3 1]        
                    for Feat=1:size(FMtmp,3)
                        FMtmp(Samp,:,Feat)=(FMtmp(Samp,:,Feat)-MeanStd(1,1))/MeanStd(1,2); % used with [2 3 1]                
                    end
                end
            else % if only 2D Matrix from FET
                FMtmp=permute(reshape(cell2mat(FMtmp).',numel(FMtmp{1}),size(FMtmp,2),[]),[2,1,3]);                
                for Feat=1:size(FMtmp,3) 
                    FMtmp(:,:,Feat)=(FMtmp(:,:,Feat)-MeanStd(Feat,1))/MeanStd(Feat,2); % used with [3 2 1]                
                end  
            end

%     % scale to  range 0-1 
%             if size(FMtmp,3)~=1 % only one feature e.g. ECG 
%                 flattenedM=permute(FMtmp,[1 3 2]);   %3D to 2D
%                 flattenedM=flattenedM(:);            %2D to 1D
%                 MInMAx(1,1)=min(flattenedM);
%                 MInMAx(1,2)=max(flattenedM);   
%             else
%                 for F=1:size(FMtmp,1) % multible features 
%                     MInMAx(F,1)=min(FMtmp(F,:));
%                     MInMAx(F,2)=max(FMtmp(F,:));
%                 end                
%             end
%             
%             if strcmp(Datapack(DT,:),'FET')~=1 % all else than FET
%                 for Samp=1:size(FMtmp,1) % used with [2 3 1]        
%                     for Feat=1:size(FMtmp,3)
%                         FMtmp(Samp,:,Feat)=(FMtmp(Samp,:,Feat)-MInMAx(1,1))./(MInMAx(1,2)-MInMAX(1,1)); % used with [2 3 1]                
%                     end
%                 end
%             else % if only 2D Matrix from FET
%                 for Feat=1:size(FMtmp,3) 
%                     FMtmp(:,:,Feat)=(FMtmp(:,:,Feat)-MInMAx(1,1))./(MInMAx(1,2)-MInMAX(1,1)); % used with [3 2 1]                
%                 end  
%             end
%             
            FeatureMatrix=FMtmp;
            

            if saving
                Saving_Session_F(FeatureMatrix,savefolderSession, Neonate,i)
                Saving_Session_A(Annotations,savefolderSession, Neonate,i)
            end

            FeatureMatrix_tmp=[FeatureMatrix_tmp ; FeatureMatrix];
            Annottmp=[Annottmp Annotations]; % 

        end % Session       
        Annotations=Annottmp;
        FeatureMatrix=FeatureMatrix_tmp;

        if saving
            Saving_F(FeatureMatrix,savefolder, Neonate,win)
            Saving_A(Annotations,savefolder, Neonate,win)
        end
        FeatureMatrix=[];
        Annotations=[];
    end %Session
end % DAtapack


%% Nested saving
    function Saving_F(FeatureMatrix,savefolder, Neonate,win)
        if exist('FeatureMatrix','var')==1
            name=inputname(1); % variable name of function input
            save([savefolder name '_' num2str(Neonate) '_win_' num2str(win)],'FeatureMatrix')
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
    
    function Saving_A(Annotations,savefolder, Neonate,win)
        if exist('Annotations','var')==1
            name=inputname(1); % variable name of function input
            save([savefolder name '_' num2str(Neonate) '_win_' num2str(win)],'Annotations')
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