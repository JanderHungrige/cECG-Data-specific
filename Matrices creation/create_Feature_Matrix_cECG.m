% function create_Feature_Matrix_cECG(Neonate)
% In the first part we add first all Time domain features to
% Feature_Matrix.
% Each Row is one Feature. Where each session of the same Feature
% is added to that row.
% In the second part we add each feature per session to each other. Each
% Feature gets it own row in the Feature_Matrix

clear
clc

Pat=[4,5,6,7,9,10,11,12,13];
path='E:\';
datapath=[path 'cECG_study\C_Processed_Data\HRV_features\'];
TFeature_path=[datapath 'timedomain\'];
FFeature_path=[datapath 'freqdomain\'];
NLFeature_path=[datapath 'nonlinear\'];

savefolder= ([path 'cECG_study\C_Processed_Data\Matrices\']);
saving=1

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
    'SEAUC'...
    };

for N=1:length(Pat)
    Neonate=Pat(N);
    FeatureMatrix=[];
%--------------- Per Patient ---------------     
    
    % all from on patient MIXED sessions TIME DOMAIN
    for j=1:length(Featurenames_time) 
        Dateien=dir([TFeature_path Featurenames_time{j,1} '_Session_*_' num2str(Neonate) '.mat']); 
        tmp=[];
        for i=1:length(Dateien)
            Dateien(i).name
            load([TFeature_path Dateien(i).name]);
            if iscell(Feature)==1 % Some features are saved in cell some as double
            tmp= [tmp, cell2mat(Feature)]; % adding each session of one feature to one single line of Feature 
            else
            tmp= [tmp, Feature]; % adding each session of one feature to one single line of Feature 
            end           
        end
        FeatureMatrix=[FeatureMatrix ;tmp]; % Adding the long single line of Feature into the Feature Matrix, where each row is one Feature

    end

    % all from on patient MIXED sessions FREQUENCY DOMAIN
    for j=1:length(Featurenames_frequency) 
        Dateien=dir([FFeature_path Featurenames_frequency{j,1} '_Session_*_' num2str(Neonate) '.mat']); 
        tmp=[];
        for i=1:length(Dateien) % adding each individual Feature
            Dateien(i).name
            load([FFeature_path Dateien(i).name]);
            if iscell(Feature)==1 % Some features are saved in cell some as double        
                tmp=[tmp , cell2mat(Feature)];% adding each session of one feature to one single line of Feature 
            else
                tmp= [tmp, Feature]; % adding each session of one feature to one single line of Feature         
            end
        end
        FeatureMatrix=[FeatureMatrix ;tmp]; % Adding the long single line of Feature into the Feature Matrix, where each row is one Feature
    end
    
    % all from on patient MIXED sessions NON LINEAR
    for j=1:length(Featurenames_nonlinear) 
        Dateien=dir([NLFeature_path Featurenames_nonlinear{j,1} '_Session_*_' num2str(Neonate) '.mat']); 
        tmp=[];
        for i=1:length(Dateien)
            Dateien(i).name
            load([NLFeature_path Dateien(i).name]);
            if iscell(Feature)==1 % Some features are saved in cell some as double
            tmp= [tmp, cell2mat(Feature)]; % adding each session of one feature to one single line of Feature 
            else
            tmp= [tmp, Feature]; % adding each session of one feature to one single line of Feature 
            end           
        end
        FeatureMatrix=[FeatureMatrix ;tmp]; % Adding the long single line of Feature into the Feature Matrix, where each row is one Feature

    end    
    
    if saving
        Saving(FeatureMatrix,savefolder, Neonate)
    end
    FeatureMatrix=[];
%--------------- Per Patient & Session --------------- 

%FOR EACH SESSION. PUT THE SESSIONS IN DIFFERENT FOLDERS AND DO ALMOST THE
%SAME AS ABOVE. SHOULD BE AROUND 4 FOLDERS

end


%% Nested saving
    function Saving(FeatureMatrix,savefolder, Neonate)
        if exist('FeatureMatrix','var')==1
            name=inputname(1); % variable name of function input
            save([savefolder name '_' num2str(Neonate) '_win_30'],'FeatureMatrix')
        else
            disp(['saving of ' name ' not possible'])
        end       
    end
    
    function Saving_Session(FeatureMatrix,savefolder, Neonate,Session)
        if exist('FeatureMatrix','var')==1
            name=inputname(1); % variable name of function input
            save([savefolder name '_' num2str(Neonate) '_Session_' num2str(Session) ],'FeatureMatrix')
        else
            disp(['saving of ' name ' not possible'])
        end       
    end
    



% end