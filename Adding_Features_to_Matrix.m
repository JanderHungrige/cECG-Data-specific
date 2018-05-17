% adding a feature to the Matrix
%*********************************************
%BACKUP YOUR FEATURE MATRIX BEFORE USING THIS
%*********************************************
clear 
clc

SignalType='cECG'; % ECG or cECG
FeatureType='nonlinear' ; % timedomain, freqdomain nonlinear
Featurename='LZECG'; % name of the feature you wnat to add to the featureMatrix

pat=[4,5,6,7,9,10,11,12,13]; 

% pat=5

if strcmp(SignalType,'ECG')
    MatrixPath= ('E:\cECG_study\C_Processed_Data\Matrices\Sessions');
    FeaturePath=(['E:\cECG_study\C_Processed_Data\HRV_features\' FeatureType]);
elseif strcmp(SignalType,'cECG')  
    MatrixPath= ('E:\cECG_study\C_Processed_Data\cMatrices\Sessions');
    FeaturePath=(['E:\cECG_study\C_Processed_Data\cHRV_features\' FeatureType]);    
end
savefolder=[MatrixPath '\'] ;
saving=1;


for i=1:length(pat)
    disp(['** working on patient ' num2str(pat(i)) '**'])
    Sessions=dir([MatrixPath '\FeatureMatrix_' num2str(pat(i)) '_*']);
    FeatureFile=dir([FeaturePath '\' Featurename '_Session_*_' num2str(pat(i)) '.mat']);
    for j=1:length(Sessions)
        disp(['Working on Session ' num2str(j)])
        load([MatrixPath '\FeatureMatrix_' num2str(pat(i)) '_Session_' num2str(j) '.mat'])
        if isempty(FeatureMatrix) % If there are many nans in the Features, it can happen, that the whole MAtrix for a session is empty. Then Skip it
            disp('FMatrix  empty')  
            continue
        end
        disp(['FMatrix size: ' num2str(length(FeatureMatrix(:,1)))])
        
        load([FeaturePath '\' FeatureFile(j).name])
        if iscolumn(Feature)
            Feature=Feature';
        end
        
        if iscell(Feature) % Some features are saved in cell some as double
            Feature=cell2mat(Feature); % adding each session of one feature to one single line of Feature 
        end
        
        [r,c]=find(any(isnan(Feature),1));
        Feature(:,c)=[];
        
        if length(FeatureMatrix(1,:))<length(Feature)
            Feature=Feature(1,1:length(FeatureMatrix(1,:)));
        end
        
        FeatureMatrix=[FeatureMatrix; Feature];
        if saving
            Saving_Session_F(FeatureMatrix,savefolder,pat(i),j)
        end
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
    