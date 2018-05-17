% This m file generates a annotation Matrix (1D Array) per patient. The
% files are merged together from each session. 
clear
clc
tic

pat=[4,5,6,7,9,10,11,12,13]; 
% pat=7
saving=1;
win=30;

%      '	1=	ActiveSleep';...
%         '	2=	QuietSleep';...
%         '	3=	Wake';...
%         '	4=	CareTaking';...
%         '	5=	UnknownBedState'...
%           6=  Transition

path='E'; % the HDD file with the patient data. Needed for loading ECG
loadfolder= [path ':\cECG_study\C_Processed_Data\Annotations\'];
savefolder= [path ':\cECG_study\C_Processed_Data\Matrices\'];


 for I=1:length(pat)
    disp('***************************************')
    disp(['Working on patient ' num2str(pat(I))])
    Neonate=pat(I);      
    Annotations=[];
  
%% ************ Load data **************
    dateien=dir([loadfolder 'Annotations_*_win_'  num2str(win) '_*_' num2str(pat(I)) '.mat']);
    for K=1:length(dateien)
        tmp=load([loadfolder dateien(K).name]);
        Annotations=[Annotations tmp.Feature];  
    end
    
     if saving
         Saving(Annotations,savefolder, Neonate,win)
     end   

 end
 
     function Saving(Annotation,savefolder, Neonate, win)
                 name=inputname(1); % variable name of function input
        if exist('Annotation','var')==1
            save([savefolder name '_win_' num2str(win) '_' num2str(Neonate)],'Annotation')
        else
            disp(['saving of ' name ' not possible'])
        end       
    end