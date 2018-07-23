%Calling RHV analysis functions

% About Patient Nr4. Session 2 (1341399361) does not have Intellivue data. Therefore,
% create DAQ data for that particular session(or pat 4 in total) rename it
% to Intellivue manually and do the same with the annotations(in total 4
% delete the others). Then you can create the matrix without lost data(6h).
% You cannot simply use the DAQ data as there are annotations missing and
% to correct for that is more difficult. 
%
% For single addition: outcomment the for loop command and end command for the sessions loop ( Line 75) and fill in e.g. S=2

clear
clc
tic


dataset='ECG'; % cECG or ECG later maybe MMC and InnerSence

filetype='DAQ'; % PAtietn 12 misses one annotation file
filetype='Intellivue'; % Patient 4 misses one Intellivue File
    
pat=[4,5,6,7,9,10,11,12,13]; 
pat=[7,9,10,11,12,13];
pat=[4,5,6]
pat=[4,5,6,7,11,13];
pat=[12]; 
pat=[4,5,6,7,9,10,11,12,13]; 

RRMethod='R'; %M or R to calculate the RR with Michiel or Ralphs algorythm
Annotators='B3A';% As Kappa is 1 we can use either of the annotations.
saving=1;
plotting=0;
win=300;
faktor=30; % how much is the data moving forward? 30s is classic
faktor_annot=1; % If we choose the Annotations per second use 30. If you use the annotations already in 30s, use 1
FS_ecg=500;


Matlabbase='C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\HRV feature creation\';

addpath(Matlabbase)
addpath('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific')
addpath('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\R peak detection')
addpath('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Annotation')
addpath('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\ECG feature creation')


loadfolder='\\code1\storage\2012-0194_neonatal_data\cECG study\RAW DATA\For_Quick_Annotator\';
% annotationfolder='\\code1\storage\2012-0194_neonatal_data\cECG study\Annotations\Data and annotations used by Jan (Bea)\Annotation\participant';
annotationfolder='C:\Users\310122653\Documents\PhD\Article_3_(cECG)\Raw Data\Annotation\participant';
savefolder=('C:\Users\310122653\Documents\PhD\Article_3_(cECG)\Processed Data\');

SavefolderAnnotations=([savefolder 'Annotations\']);


 for I=1:length(pat)
    disp('***************************************')
    disp(['Working on patient ' num2str(pat(I))])
    Neonate=pat(I);      
    
  
%% ************ Load data **************
    
    if strcmp('cECG',dataset)==1 % as the cECG values are stored with the DAQ
        filetype='DAQ'; % PAtietn 12 misses one annotation file
    end

%     filetype='DAQ';nn% use only when Pat=4 and S=2    

    filelocation=([loadfolder 'participant' num2str(Neonate) '\' filetype '\']);
    Datum=dir([filelocation '2012*']); % folder name differnt on this stage
    filelocation=([filelocation Datum.name '\']);
    Sessions=dir([filelocation filetype '*']);
    
    for S=1:length(Sessions)

%% ************ Load annotations (1s) **************  
%loading 1 secondannotations for this particular patient/session

        Annotation=loading_annotations(Neonate,Sessions(S,1).name,annotationfolder);
        disp('* Annotation loaded')

 %% ************ Window  ECG /  Annotation signals 
      
        [Annotations_win_300, Annotations_win_30]=SlidingWindow_Annotations(Annotation,Neonate,saving,SavefolderAnnotations,win,S,faktor_annot);           
        
        
    end %Sessionp

 end% Patient
toc
%% Nested saving
    function Saving(Feature,savefolder, Neonate, win,S)
        if exist('Feature','var')==1
            name=inputname(1); % variable name of function input
            save([savefolder name '_Session_' num2str(S) '_pat_' num2str(Neonate)],'Feature')
        else
            disp(['saving of ' name ' not possible'])
        end       
    end
 