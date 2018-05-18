%Calling RHV analysis functions

% About Patient Nr4. Session 2 (1341399361) does not have Intellivue data. Therefore,
% create DAQ data for that particular session(or pat 4 in total) rename it
% to Intellivue manually and do the same with the annotations(if total 4
% delete the others). Then you can create the matrix without lost data(6h).
% You cannot simply use the DAQ data as there are annotations missing and
% to correct for that is more difficult. 
%
% For single addition: outcomment the for loop command and end command for the sessions loop ( Line 75) and fill in e.g. S=2

clear
clc
tic


dataset='ECG'; % cECG or ECG later maybe MMC and InnerSence

pat=[4,5,6,7,9,10,11,12,13]; 
pat=[7,9,10,11,12,13];
pat=[4,5,6]
pat=[4,5,6,7,11,13];


RRMethod='R'; %M or R to calculate the RR with Michiel or Ralphs algorythm
saving=1;
plotting=0;
win=300;
faktor=30; % how much is the data moving forward? 30s is classic
FS_ecg=500;

Matlabbase='C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\HRV feature creation\';

addpath(Matlabbase)
addpath('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific')
addpath('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\R peak detection')
addpath('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\Annotation')
addpath('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\ECG feature creation')

path='E'; % the HDD file with the patient data. Needed for loading ECG
loadfolder=([path ':\cECG_study\B_Annotations\Datafiles\For Quick Annotator\']);
% savefolder= ([path ':\cECG_study\C_Processed_Data\']);
savefolder=('C:\Users\310122653\Documents\PhD\Article_3_(cECG)\Processed Data\');
% savefolder=('E:\cECG_study\C_Processed_Data\For Xi\');

SavefolderAnnotations=([ savefolder 'Annotations\']);

if strcmp(RRMethod,'R')
    if strcmp('ECG',dataset)==1
        savefolderHRVtime= ([ savefolder 'HRV_features\timedomain\']);
        savefolderHRVfreq= ([ savefolder 'HRV_features\freqdomain\']);        
        savefolderHRVnonlin= ([ savefolder 'HRV_features\nonlinear\']);
        savefolderECG= ([ savefolder 'HRV_features\ECG\']);
        savefolderEDR=([ savefolder 'HRV_features\EDR\']);
        savefolderRR=([ savefolder 'HRV_features\RR\']);
        savefolderResp=([ savefolder 'HRV_features\Resp\']);
        savefolderAnnotations=([ savefolder 'HRV_features\Annotations\']);

    elseif strcmp('cECG',dataset)==1
        savefolderHRVtime= ([ savefolder 'cHRV_features\timedomain\']);
        savefolderHRVfreq= ([ savefolder 'cHRV_features\freqdomain\']);        
        savefolderHRVnonlin= ([ savefolder 'cHRV_features\nonlinear\']);
        savefolderECG= ([ savefolder 'cHRV_features\ECG\']);
        savefolderEDR=([ savefolder 'cHRV_features\EDR\']);
        savefolderRR=([ savefolder 'cHRV_features\RR\']);
        savefolderResp=([ savefolder 'cHRV_features\Resp\']); 
        savefolderAnnotations=([ savefolder 'cHRV_features\Annotations\']);
        
    else
        disp('Error: wrong dataset string. Line 7 in CallingHRVfunctions_for_cECG')
        stop
    end 
elseif  strcmp(RRMethod,'M')
    if strcmp('ECG',dataset)==1
        savefolderHRVtime= ([ savefolder 'HRV_features\timedomainM\']);
        savefolderHRVfreq= ([ savefolder 'HRV_features\freqdomainM\']);        
        savefolderHRVnonlin= ([ savefolder 'HRV_features\nonlinearM\']);
        savefolderECG= ([ savefolder 'HRV_features\ECGM\']);
        savefolderEDR=([ savefolder 'HRV_features\EDRM\']);
        savefolderRR=([ savefolder 'HRV_features\RRM\']);
        savefolderResp=([ savefolder 'HRV_features\RespM\']);
        
    elseif strcmp('cECG',dataset)==1
        savefolderHRVtime= ([ savefolder 'cHRV_features\timedomainM\']);
        savefolderHRVfreq= ([ savefolder 'cHRV_features\freqdomainM\']);        
        savefolderHRVnonlin= ([ savefolder 'cHRV_features\nonlinearM\']); 
        savefolderECG= ([ savefolder 'cHRV_features\ECGM\']);
        savefolderEDR=([ savefolder 'cHRV_features\EDRM\']);
        savefolderRR=([ savefolder 'cHRV_features\RRM\']);
        savefolderResp=([ savefolder 'cHRV_features\RespM\']);
        
    else
        disp('Error: wrong dataset string. Line 7 in CallingHRVfunctions_for_cECG')
        stop
    end 
end


 for I=1:length(pat)
    disp('***************************************')
    disp(['Working on patient ' num2str(pat(I))])
    Neonate=pat(I);      
    
  
%% ************ Load data **************

%     filetype='DAQ'; % PAtietn 12 misses one annotation file
    filetype='Intellivue'; % Patient 4 misses one Intellivue File
    
    if strcmp('cECG',dataset)==1 % as the cECG values are stored with the DAQ
        filetype='DAQ'; % PAtietn 12 misses one annotation file
    end

%     filetype='DAQ';nn% use only when Pat=4 and S=2    

    filelocation=([loadfolder 'participant' num2str(Neonate) '\' filetype '\']);
    Datum=dir([filelocation '2012*']); % folder name differnt on this stage
    filelocation=([filelocation Datum.name '\']);
    Sessions=dir([filelocation filetype '*']);
    
    for S=1:length(Sessions)
%         S=8
%         if pat(I)==4 && S>1 %For patient 4 the Second Session has to be taken from DAQ. To not rename all manually shift over S=2
%             S=S+1; %Leaving Session 2 empty. That has to be added "manually" using the DAQ. Remove the for loop and this if as well as the if before the end. Then set S=2
%         end
%         if pat(I)==4 && S==2
%             continue
%             %For patient 4 the Second Session has to be taken from DAQ insead of Intellivue as the Intellivue file is missing. To not rename all manually shift over S=2
%             %Leaving Session 2 empty. That has to be added "manually" using the DAQ. Put filetype='DAQ', Remove this if and set S=2
%         end
%                 

        disp('- - - - - - - - - - - - - - - - - - - ')
        disp(['Working on session: ' Sessions(S,1).name ' NR:' num2str(S) '/' num2str(length(Sessions))])
        load([filelocation Sessions(S,1).name])


%% ************ Load annotations (1s) **************  
%loading 1 secondannotations for this particular patient/session

        Annotation=loading_annotations(Neonate,Sessions(S,1).name);
        disp('* Annotation loaded')

 %% ************ Window  ECG /  Annotation signals 

        if strcmp('cECG',dataset)==1
            t_ECG=linspace(0,length(cECG.values)/FS_ecg,length(cECG.values))';
            [ECG_win_300,ECG_win_30,t_ECG_300,t_ECG_30]=SlidingWindow_ECG(cECG.values,t_ECG,Neonate,saving,folder,faktor,win,S); 
            [Annotations_win_300, Annotations_win_30]=SlidingWindow_Annotations(Annotation,t_ECG,Neonate,saving,savefolder,win,S,faktor);
        elseif strcmp('ECG',dataset)==1
            t_ECG=linspace(0,length(ECG.values)/FS_ecg,length(ECG.values))';
            t_Resp=linspace(0,length(Resp.values)/FS_ecg,length(Resp.values))';
            t_EDR=linspace(0,length(EDR.values)/FS_ecg,length(EDR.values))';            
            % The differnec in t_300 and t_ECG_300 is that t_ECG_300 is a
            % continuous run of time, while t_300 is 0 to t for each cell element
           [ECG_win_300,ECG_win_30,t_ECG_300,t_ECG_30]=SlidingWindow_ECG(ECG.values,t_ECG,Neonate,saving,savefolderEDR,faktor,win,S); 
           [Annotations_win_300, Annotations_win_30]=SlidingWindow_Annotations(Annotation,t_ECG,Neonate,saving,savefolderAnnotations,win,S,faktor);           
%            [Resp_win_300,Resp_win_30,t_Resp_300,t_Resp_30]=SlidingWindow_Resp(Resp.values,t_Resp,Neonate,win,saving,savefolderResp,Sessions(S,1).name,faktor,S); 
           [EDR_win_300,EDR_win_30,t_EDR_300,t_EDR_30]=SlidingWindow_EDR(EDR.values,t_EDR,Neonate,saving,savefolderEDR,faktor,win,S); 
             disp(['* Data is merged into windows of length: ' num2str(win) 's and ' num2str(30) 's'] )  
        else
            disp('Error. no cECG or ECG could be found: Calling RHV analysis functions; line 75-86')
        end
    %% ************ Creating RR signal for ECG-Signal **************
        Ralphsfactor={1;-1;-1;1;1;-1;1;-1; 1; 1;-1; 1;-1; 1;-1;-1;-1;-1};%Determine if the ECG signal should be turned -1 or not 1. 
                     %1  2 3  4 5  6 7  8  9  10 11 12 13 14 15 16 17 18
        padding=0; %Determine if the RR should be same length as ECG. Don`t have to be
        plotting=0; %plotting Ralphs RR detection
        
        if strcmp(RRMethod,'M')
%Michiel
%             for R=1:length(ECG_win_300)
%                 t_300{1,R}=linspace(0,length(ECG_win_300{1,R})/FS_ecg,length(ECG_win_300{1,R}))';
%                 if all(isnan(ECG_win_300{1,R}'))==1 || range(ECG_win_300{1,R}')==0  % if only Nan Ralph cannot handle it or if all values are the same (Flat line)
%                      RR_300{R,1}=NaN(1,length(ECG_win_300{1,R})) ;
%                 else
%                     ECG_win_300{1,R}(isnan(ECG_win_300{1,R})) = []; %Michiels cannot handle a single NAN
%                     t_300{1,R}=linspace(0,length(ECG_win_300{1,R})/FS_ecg,length(ECG_win_300{1,R}))';                   
%                     [RR_idxM] = streamingpeakdetection(ECG_win_300{1,R}', 500, [60 256], plotting, 18.5, 1024);
%                     RR_300{R,1}=diff(RR_idxM.peakPositionArray)./500; % Calculating the time between the R peaks in seconds
%                     RR_300{R,1}=[NaN RR_300{R,1}];% NAn is needed to make RR and RR_idx same length which is needed for lomb scargle
%                     RR_idx_300{R,1}=RR_idxM.peakPositionArray; 
%                 end
%             end           
%             for R=1:length(ECG_win_30)
%                 t_30{1,R}=linspace(0,length(ECG_win_30{1,R})/FS_ecg,length(ECG_win_30{1,R}))';
%                 if all(isnan(ECG_win_30{1,R}'))==1 || range(ECG_win_30{1,R}')==0  % if only Nan Ralph cannot handle it or if all values are the same (Flat line)
%                      RR_30{R,1}=NaN(1,length(ECG_win_30{1,R})) ;
%                 else
%                     ECG_win_30{1,R}(isnan(ECG_win_30{1,R})) = []; %Michiels cannot handle a single NAN
%                     t_30{1,R}=linspace(0,length(ECG_win_30{1,R})/FS_ecg,length(ECG_win_30{1,R}))';                   
%                     [RR_idxM] = streamingpeakdetection(ECG_win_30{1,R}', 500, [60 256], plotting, 18.5, 1024);
%                     RR_30{R,1}=diff(RR_idxM.peakPositionArray)./500; % Calculating the time between the R peaks in seconds
%                     RR_30{R,1}=[NaN RR_30{R,1}];                    
%                     RR_idx_30{R,1}=RR_idxM.peakPositionArray; 
%                 end
%             end                
%         elseif strcmp(RRMethod,'R')
% %Ralph            
%             for R=1:length(ECG_win_300)
%                 t_300{1,R}=linspace(0,length(ECG_win_300{1,R})/FS_ecg,length(ECG_win_300{1,R}))';
%                 if all(isnan(ECG_win_300{1,R}))==1 || range(ECG_win_300{1,R})==0  % if only Nan Ralph cannot handle it or if all values are the same (Flat line)
%                    RR_300{R,1}=NaN(1,length(ECG_win_300{1,R})) ;
%                 else
%                     [RR_idx_300{R,1}, ~, ~, ~, ~, RR_300{R,1}, ~] = ecg_find_rpeaks(t_300{1,R},Ralphsfactor{Neonate,1}*ECG_win_300{1,R}, FS_ecg, 250,plotting,0); %, , , maxrate,plotting,saving   -1* because Ralph optimized for a step s slope, we also have steep Q slope. inverting fixes that probel 
%                 end
%             end
%             for R=1:length(ECG_win_30)  
%                 t_30{1,R}=linspace(0,length(ECG_win_30{1,R})/FS_ecg,length(ECG_win_30{1,R}))';        
%                 if all(isnan(ECG_win_30{1,R}))==1 || range(ECG_win_30{1,R})==0 % if all elements are NAN or the same value, R peaks cannot be calculated
%                    RR_30{R,1}=NaN(1,length(ECG_win_30{1,R})) ;
%                 else        
%                 [RR_idx_30{R,1}, ~, ~, ~, ~, RR_30{R,1}, ~] = ecg_find_rpeaks(t_30{1,R},Ralphsfactor{Neonate,1}*ECG_win_30{1,R}, FS_ecg, 250,plotting,0); %, , , maxrate,plotting,saving   -1* because Ralph optimized for a step s slope, we also have steep Q slope. inverting fixes that probel             
%                 end
%             end        
%         end
%         disp('* RR calcuated')
%         if saving; RR=RR_30 ;  Saving(RR,savefolderRR, Neonate, win,S);end
    %% ************ Creating spectrum for ECG-Signal **************         
%        [powerspectrum,f]=Lomb_scargel_single(RR_300,RR_idx_300,t_300,Neonate,saving,savefolderHRVfreq,win,S) ;
%         disp('* Periodogram calculated')
% 
%         
%     %% ************ CALCULATE FEATURES **************
% 
%     %%%%%%%% ECG TIME DOMAIN     
%         disp('ECG time domain analysis start') 
% 
%         Beats_per_Epoch(RR_300,Neonate,saving,savefolderHRVtime,win,S)   % S for session number
%              disp('- BpE finished')
%         linelength(ECG_win_300,t_300,Neonate,saving,savefolderHRVtime,win,S)     
%              disp('- Linelength finished')
%         meanarclength(ECG_win_30,t_30,Neonate,saving,savefolderHRVtime,win,faktor,S) 
%              disp('- Mean linelength finished')
%         SDLL(ECG_win_30,t_30,Neonate,saving,savefolderHRVtime,win,faktor,S) %Standart derivation of 5min linelength
%              disp('- SDLL finsihed')
%         SDaLL(ECG_win_30,t_30,Neonate,saving,savefolderHRVtime,win,faktor,S) %Standart derivation of 30s linelength meaned over 5min
%              disp('- SDaLL finished') 
% 
% 
%   %%%%%%%% HRV TIME DOMAIN
%         disp('HRV time domain analysis start')
% 
%         SDNN(RR_300,Neonate,saving,savefolderHRVtime,win,S);
%             disp('- SDNN finished') 
%         RMSSD(RR_300,Neonate,saving,savefolderHRVtime,win,S);
%             disp('- RMSSD finished')  
%         NNx(RR_300,Neonate,saving,savefolderHRVtime,win,S);
%             disp('- NNx finished') 
%         pNNx(RR_300,Neonate,saving,savefolderHRVtime,win,S);
%             disp('- pNNx finished') 
%         SDANN(RR_30,Neonate,saving,savefolderHRVtime,win,faktor,S);
%             disp('- SDANN finished')
%         pDec(RR_300,Neonate,saving,savefolderHRVtime,win,S);
%             disp('- pDEC finished') 
%         SDDec(RR_300,Neonate,saving,savefolderHRVtime,win,S);
%            disp('- SDDec finished')
%  
%           
%   %%%%%%% HRV Frequency domain
%         disp('Frequency time domain start')
% 
%         freqdomainHRV (powerspectrum,f,Neonate,win,saving,savefolderHRVfreq,S)
%            disp('- Frequency finished') 
% 
% 
%     %%%%%% HRV Non linear
%         disp('Nonlinear analysis start')
% 
%         SampEn_QSE_SEAUC(RR_300,Neonate,saving,savefolderHRVnonlin,win,faktor,S ) %
%             disp('- SampEn QSE SEAUCfinished')
%         LempelZivECG(ECG_win_300,Neonate,saving,savefolderHRVnonlin,win,S)  
%           disp('- LepelZiv ECG finished')         
%         LempelZivRR(RR_300,Neonate,saving,savefolderHRVnonlin,win,S);
%           disp('- LepelZiv HRV finished')   


        clearvars ECG_win_300 ECG_win_30 t_ECG_300 t_ECG_30 RR_idx_300 RR_300 RR_idx_30 RR_30 powerspectrum f
        
    end %Sessionp

 end% Patient
toc
%% Nested saving
    function Saving(Feature,savefolder, Neonate, win,S)
        if exist('Feature','var')==1
            name=inputname(1); % variable name of function input
            save([savefolder name '_Session_' num2str(S) 'pat_' num2str(Neonate)],'Feature')
        else
            disp(['saving of ' name ' not possible'])
        end       
    end
 