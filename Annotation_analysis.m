
clear
clc
tic


%      '	1=	ActiveSleep';...
%         '	2=	QuietSleep';...
%         '	3=	Wake';...
%         '	4=	CareTaking';...
%         '	5=	UnknownBedState'...
%           6=  Transition

Matlabbase='C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\HRV feature creation\';

addpath(Matlabbase)
addpath('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific')
addpath('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\R peak detection')
addpath('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Annotation')
addpath('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\ECG feature creation')


path='E'; % the HDD file with the patient data. Needed for loading ECG
loadfolder=([path ':\cECG_study\B_Annotations\Datafiles\For Quick Annotator\']);
savefolder= ([path ':\cECG_study\C_Processed_Data\']);
savefolderHRVtime= ([path ':\cECG_study\C_Processed_Data\HRV_features\timedomain\']);
savefolderHRVfreq= ([path ':\cECG_study\C_Processed_Data\HRV-features\freqdomain\']);        



pat=[4,5,6,7,9,10,11,12,13]; 
% pat=7
saving=0;
plotting=0;
win=300;
winXi=300;
faktor=30; % how much is the data moving forward? 30s is classic
FS_ecg=500;

totalannot=[]

 for I=1:length(pat)
    disp('***************************************')
    disp(['Working on patient ' num2str(pat(I))])
    Neonate=pat(I);      
    
  
%% ************ Load data **************

    filetype='Intellivue';
    filelocation=([loadfolder 'participant' num2str(Neonate) '\' filetype '\']);
    Datum=dir([filelocation '2012*']); % folder name differnt on this stage
    filelocation=([filelocation Datum.name '\']);
    Sessions=dir([filelocation filetype '*']);
    
    for S=1:length(Sessions)
        disp('- - - - - - - - - - - - - - - - - - - ')
        disp(['Working on session: ' Sessions(S,1).name])
        load([filelocation Sessions(S,1).name])
        
        
%% ************ Load annotations (1s) **************  
%loading 1 secondannotations for this particular patient/session

    Annotation=loading_annotations(Neonate,Sessions(S,1).name);
    totalannot=[totalannot ; Annotation]; % collecting all annotations for analysis
    end    

 end
 
    [counts,centers]=hist(totalannot,6);
    figure, bar([1,2,3,4,5,6],counts)
    
    countspercent=100/length(totalannot)*counts ;
    figure, bar([1,2,3,4,5,6],countspercent)
    set(gcf,'color','w'); %background white
    l{1}='AS';l{2}= 'QS'; l{3}= 'Wake'; l{4}= 'CareT' ; l{5}= 'NA'; l{6}= 'Trans';
    set(gca,'xticklabel', l) % CHange X tick names
    labels = arrayfun(@(countspercent) num2str(countspercent ,'%2.1f'),countspercent,'UniformOutput',false); % get value in graph
    text([1,2,3,4,5,6],countspercent,labels ,...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom') 
    