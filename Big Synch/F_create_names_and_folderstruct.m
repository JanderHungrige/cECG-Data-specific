% Create the names and folder stucture for the Quick annotator tool


clc 
clear
clc
%% Declarations
% pat=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];

cECG=1; % are there cECG files? 
video=1; %If you want to create also video = 1 otherwise (already existing/time saving) then = 0

if cECG==1
    pat=[1,2,3,5,7,8,9,10,11,12,13,14]; %cECG
elseif cECG==0
    pat=[4,6,14,15,16,17,18];%no cECG
end
%% 
pat=[9];
missing_synched_data=[];

Matlabfolder=('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific');
QAfolder=([Matlabfolder '\Quick Annotator']);
path='e:\';

folder=dir([path 'cECG_study\A_RawData\*_test*' ]);



for i=1:length(pat)
  disp(['working on patient ' num2str(pat(1,i))])
        Neonate=pat(i);
        cd(QAfolder)
    if cECG==1
        create_folderandName_structur(Neonate,video)
    elseif cECG==0
        create_folderandName_structur_no_cECG(Neonate,path,QAfolder,video)   
    end
    
end % for each patient


%  missing_synched_data=missing_synched_data(~cellfun('isempty',missing_synched_data))  ;
%  missing_cECG_data=missing_cECG_data(~cellfun('isempty',missing_cECG_data))  ;
%  
% disp('Finished')
% disp(' synch folders are missing for: ')
% % for j=1:length(missing_synched_data)
%     disp(missing_synched_data)
% disp(' cECG folders are missing for: ')
%     disp(missing_cECG_data)
%     