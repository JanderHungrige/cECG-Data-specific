
%Calling RHV analysis functions

% About Patient Nr4. Session 2 (1341399361) does not have Intellivue data. Therefore,
% create DAQ data for that particular session(or pat 4 in total) rename it
% to Intellivue manually and do the same with the annotations(if total 4
% delete the others). Then you can create the matrix without lost data(6h).
% You cannot simply use the DAQ data as there are annotations missing and
% to correct for that is more difficult. 

%   THIS IS DONE BY JUST COPYING THE RAW 2 SESSION FROM daq INTO THE RAW
%   DATA FOLDER OF INTELLIVUE (ALREADY DONE; Check nevertheless)
%
% For single addition: outcomment the for loop command and end command for the sessions loop ( Line 75) and fill in e.g. S=2

clear
clc
tic

shutItDown=1; % if you want to shut the PC down after finishing all calculations;

dataset='ECG'; % cECG or ECG later maybe MMC and InnerSence
filetype='DAQ'; % PAtietn 12 misses one annotation file
filetype='Intellivue'; % Patient 4 misses one Intellivue File
    
pat=[4,5,6,7,9,10,11,12,13]; 
% pat=[4,5,6,7,11,13];
% pat=[7,11]
pat=[4]; 
shutdown=0
saving=1
user='c3po';

RRMethod='R'; %M or R to calculate the RR with Michiel or Ralphs algorythm
Annotators='B3A';% As Kappa is 1 we can use either of the annotations.
plotting=0;
win=300;
faktor=30; % how much is the data moving forward? 30s is classic
faktor_annot=1; % If we choose the Annotations per second use 30. If you use the annotations already in 30s, use 1
FS_ecg=500;

Pat_weight=[1606, 2410, 1160, 1845, 1420, 1110, 755, 2080, 2480];
Pat_GA=[30*7+3, 39*7+6, 27*7+3, 31*7+5, 30*7+3, 29*7+2, 27*7, 33*7+2, 27*7+3];
Pat_CA=[30*7+3+6, 39*7+4+6, 27*7+3+7+4, 31*7+5+3*7+1, 30*7+3+5, 29*7+2+7+2, 27*7+2*7+1, 33*7+2+7+3, 27*7+3+8*7+3];
Pat_GACA=Pat_CA-Pat_GA;
NICU_info=[1,1,1,2,1,1,1,2,1];% NICU=1, NMCU=2
C02=[2,2,2,1,1,1,1,3,2 ]; %1=CPAP, 2=lowflow 3=no

% Order AGe and weight after 10, 25 and 75 percentiles
Pat_weight(Pat_weight<=765)=1;
Pat_weight(Pat_weight>765 & Pat_weight<=1050)=2;
Pat_weight(Pat_weight>1050 & Pat_weight<=1665)=3;
Pat_weight(Pat_weight>1665)=4;     

Pat_GA(Pat_GA<=26.9143*7)=1;
Pat_GA(Pat_GA>26.9143*7 & Pat_GA<=27.8571*7)=2;
Pat_GA(Pat_GA>27.8571*7 & Pat_GA<=30.4286*7)=3;
Pat_GA(Pat_GA>30.4286*7)=4;  

Pat_CA(Pat_CA<=30.0571*7)=1;
Pat_CA(Pat_CA>30.0571*7 & Pat_CA<=31.2857 *7)=2;
Pat_CA(Pat_CA>31.2857*7 & Pat_CA<=34.6429*7)=3;
Pat_CA(Pat_CA>34.6429*7)=4;  

Pat_GACA(Pat_GACA<=0.8571*7)=1;
Pat_GACA(Pat_GACA>0.8571*7 & Pat_GACA<=1.4286  *7)=2;
Pat_GACA(Pat_GACA>1.4286 *7 & Pat_GACA<=5.8214*7)=3;
Pat_GACA(Pat_GACA>5.8214*7)=4;  

if strcmp(user,'c3po')
    basepath='C:\Users\C3PO';
elseif strcmp(user,'Philips')
    basepath='C:\Users\310122653';
end

Matlabbase=[basepath '\Documents\GitHub\cECG-Data-specific'];
cd(Matlabbase)
addpath([basepath '\Documents\GitHub\cECG-Data-specific'])
addpath([basepath '\Documents\GitHub\cECG-Data-specific\Annotation'])
addpath([basepath '\Documents\GitHub\cECG-Data-specific\Create data'])
addpath([basepath '\Documents\GitHub\Joined_Matlab'])
addpath([basepath '\Documents\GitHub\Joined_Matlab\R peak detection'])
addpath([basepath '\Documents\GitHub\Joined_Matlab\ECG feature creation'])
addpath([basepath '\Documents\GitHub\Joined_Matlab\HRV feature creation\'])


path='E'; % the HDD file with the patient data. Needed for loading ECG
% loadfolder=([path ':\cECG_study\B_Annotations\Datafiles\For Quick Annotator\']);
% annotationfolder= [path ':\cECG_study\B_Annotations\' Annotators '\participant'];
% savefolder= ([path ':\cECG_study\C_Processed_Data\']);

%Folders after HDD crash. Now with VPN from Philips storage
if strcmp(user,'c3po')
    loadfolder='D:\PhD\Article_3_(cECG)\RAW DATA\';
    annotationfolder='D:\PhD\Article_3_(cECG)\Annotation\participant';
    savefolder='D:\PhD\Article_3_(cECG)\Processed Data\';
elseif strcmp(user,'Philips')
    % annotationfolder='\\code1\storage\2012-0194_neonatal_data\cECG study\Annotations\Data and annotations used by Jan (Bea)\Annotation\participant';
    loadfolder='\\code1\storage\2012-0194_neonatal_data\cECG study\RAW DATA\For_Quick_Annotator\';
    annotationfolder='C:\Users\310122653\Documents\PhD\Article_3_(cECG)\Raw Data\Annotation\participant';
    savefolder=('C:\Users\310122653\Documents\PhD\Article_3_(cECG)\Processed Data\');
end
    
    
SavefolderAnnotations=([ savefolder 'Annotations\']); mkdir (SavefolderAnnotations) ;

if strcmp('ECG',dataset)==1
    savefolderHRVtime= ([ savefolder 'HRV_features\timedomain\']); mkdir (savefolderHRVtime) ;
    savefolderHRVfreq= ([ savefolder 'HRV_features\freqdomain\']); mkdir (savefolderHRVfreq) ;    
    savefolderHRVnonlin= ([ savefolder 'HRV_features\nonlinear\']); mkdir (savefolderHRVnonlin) ;
    savefolderECG= ([ savefolder 'HRV_features\ECG\']);mkdir (savefolderECG) ;
    savefolderEDR=([ savefolder 'HRV_features\EDR\']);mkdir (savefolderEDR) ;
    savefolderRR=([ savefolder 'HRV_features\RR\']);mkdir (savefolderRR) ;
    savefolderResp=([ savefolder 'HRV_features\Resp\']);mkdir (savefolderResp) ;
    savefolderAGEWEight=([ savefolder 'HRV_features\AGEWEight\']);mkdir (savefolderAGEWEight) ;

elseif strcmp('cECG',dataset)==1
    savefolderHRVtime= ([ savefolder 'cHRV_features\timedomain\']);
    savefolderHRVfreq= ([ savefolder 'cHRV_features\freqdomain\']);        
    savefolderHRVnonlin= ([ savefolder 'cHRV_features\nonlinear\']);
    savefolderECG= ([ savefolder 'cHRV_features\ECG\']);
    savefolderEDR=([ savefolder 'cHRV_features\EDR\']);
    savefolderRR=([ savefolder 'cHRV_features\RR\']);
    savefolderResp=([ savefolder 'cHRV_features\Resp\']); 
    savefolderAGEWEight=([ savefolder 'cHRV_features\AGEWEight\']);mkdir (savefolderAGEWEight) ;


else
    disp('Error: wrong dataset string. Line 7 in CallingHRVfunctions_for_cECG')
    stop
end 

ExtraSession=0; % for Pat4 S=2 

 for I=1:length(pat)
    disp('***************************************')

    disp(['Working on patient ' num2str(pat(I))])
    Neonate=pat(I);      
    
  
%% ************ Load data **************
    
    if strcmp('cECG',dataset)==1 % as the cECG values are stored with the DAQ
        filetype='DAQ'; % PAtietn 12 misses one annotation file
    end
    

    filelocation=([loadfolder 'participant' num2str(Neonate) '\' filetype '\']);
    Datum=dir([filelocation '2012*']); % folder name differnt on this stage
    filelocation=([filelocation Datum.name '\']);
    Sessions=dir([filelocation filetype '*']);
    
    for S=1:length(Sessions)
% Pat4 Session2 is just copied manually in the raw folder and the results should be renamed to intellivue

        disp('- - - - - - - - - - - - - - - - - - - ')
        disp(['Working on session: ' Sessions(S,1).name ' NR:' num2str(S) '/' num2str(length(Sessions))])
        load([filelocation Sessions(S,1).name])


%% ************ Load annotations (1s) **************  
%loading 1 secondannotations for this particular patient/session

        Annotation=loading_annotations(Neonate,Sessions(S,1).name,annotationfolder);
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
           [ECG_win_300,ECG_win_30,t_ECG_300,t_ECG_30]=SlidingWindow_ECG(ECG.values,t_ECG,Neonate,saving,savefolderECG,faktor,win,S); 
           [Annotations_win_300, Annotations_win_30]=SlidingWindow_Annotations(Annotation,Neonate,saving,SavefolderAnnotations,win,S,faktor_annot);           
%            [Resp_win_300,Resp_win_30,t_Resp_300,t_Resp_30]=SlidingWindow_Resp(Resp.values,t_Resp,Neonate,win,saving,savefolderResp,Sessions(S,1).name,faktor,S); 
%            [EDR_win_300,EDR_win_30,t_EDR_300,t_EDR_30]=SlidingWindow_EDR(EDR.values,t_EDR,Neonate,saving,savefolderEDR,faktor,win,S); 
             disp(['* Data is merged into windows of length: ' num2str(win) 's and ' num2str(30) 's'] )  
        else
            disp('Error. no cECG or ECG could be found: Calling RHV analysis functions; line 75-86')
        end
    %% ************ Creating RR signal for ECG-Signal **************
        Ralphsfactor={1;-1;-1;1;1;-1;1;-1; 1; 1;-1; 1;-1; 1;-1;-1;-1;-1};%Determine if the ECG signal should be turned -1 or not 1. 
                     %1  2 3  4 5  6 7  8  9  10 11 12 13 14 15 16 17 18
        padding=0; %Determine if the RR should be same length as ECG. Don`t have to be
        plotting=0; %plotting Ralphs RR detection
        

%Ralph            
        for R=1:length(ECG_win_300)
            t_300{1,R}=linspace(0,length(ECG_win_300{1,R})/FS_ecg,length(ECG_win_300{1,R}))';
            if all(isnan(ECG_win_300{1,R}))==1 || range(ECG_win_300{1,R})==0  % if only Nan Ralph cannot handle it or if all values are the same (Flat line)
               RR_300{R,1}=NaN(1,length(ECG_win_300{1,R})) ;
            else
                [RR_idx_300{R,1}, ~, ~, ~, ~, RR_300{R,1}, ~] = ecg_find_rpeaks(t_300{1,R},Ralphsfactor{Neonate,1}*ECG_win_300{1,R}, FS_ecg, 250,plotting,0); %, , , maxrate,plotting,saving   -1* because Ralph optimized for a step s slope, we also have steep Q slope. inverting fixes that probel 
            end
        end
        for R=1:length(ECG_win_30)  
            t_30{1,R}=linspace(0,length(ECG_win_30{1,R})/FS_ecg,length(ECG_win_30{1,R}))';        
            if all(isnan(ECG_win_30{1,R}))==1 || range(ECG_win_30{1,R})==0 % if all elements are NAN or the same value, R peaks cannot be calculated
               RR_30{R,1}=NaN(1,length(ECG_win_30{1,R})) ;
            else        
            [RR_idx_30{R,1}, ~, ~, ~, ~, RR_30{R,1}, ~] = ecg_find_rpeaks(t_30{1,R},Ralphsfactor{Neonate,1}*ECG_win_30{1,R}, FS_ecg, 250,plotting,0); %, , , maxrate,plotting,saving   -1* because Ralph optimized for a step s slope, we also have steep Q slope. inverting fixes that probel             
            end
        end        
        
        disp('* RR calcuated')
        if saving; RR=RR_30 ; SavingC(RR,savefolderRR, Neonate, win,Sessions(S,1).name,S); disp('* RR saved'); end
        
        
        %% ************ Creating EDR signal from 30s epoch ECG **************
        [EDR_30] =Respiration_from_ECG(ECG_win_30,RR_idx_30,RR_30,500);        
        [EDR_300]=Respiration_from_ECG(ECG_win_300,RR_idx_300,RR_300,500);   
        disp('* EDR calculated')
        if saving
            SavingC(EDR_30,savefolderEDR, Neonate, win,Sessions(S,1).name,S)
            SavingC(EDR_300,savefolderEDR, Neonate, win,Sessions(S,1).name,S)        
            disp('* EDR saved')
        end              
    %% ************ Creating spectrum for ECG-Signal **************         
       [powerspectrum,f]=Lomb_scargel_single(RR_300,RR_idx_300,t_300) ;
       [powerspectrumEDR,fEDR]=Lomb_scargel_single(EDR_300,RR_idx_300,t_300) ;
       disp('* Periodogram calculated')
        if saving    
            SavingC(powerspectrum,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S)
            SavingC(powerspectrumEDR,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S)
            disp('* Spectrums saved')
        end


    %%  ************ AGE & Weight ************** 

    
    for k=1:length(RR)
        Birthweight{k}=Pat_weight(I);
        GA{k}=Pat_GA(I); 
        CA{k}=Pat_CA(I);
        Age_diff{k}=Pat_GACA(I);
    end
    if saving
        SavingC(Birthweight,savefolderAGEWEight, Neonate, win,Sessions(S,1).name,S)
        SavingC(GA,savefolderAGEWEight, Neonate, win,Sessions(S,1).name,S)
        SavingC(CA,savefolderAGEWEight, Neonate, win,Sessions(S,1).name,S)
        SavingC(Age_diff,savefolderAGEWEight, Neonate, win,Sessions(S,1).name,S)              
        disp('* Age and Weight saved')
    end       
    clearvars Birthweight GA CA Age_diff

    %% ************ CALCULATE FEATURES **************

    %%%%%%%% ECG TIME DOMAIN     
        disp('ECG time domain analysis start') 

        BpE=Beats_per_Epoch(RR_300);SavingC(BpE,savefolderHRVtime, Neonate, win,Sessions(S,1).name,S)   % S for session number
             disp('- BpE finished')
        LL=linelength(ECG_win_300,t_300);SavingC(LL,savefolderHRVtime, Neonate, win,Sessions(S,1).name,S)   
             disp('- Linelength finished')
        aLL=meanarclength(ECG_win_30,t_30,faktor,win);SavingC(aLL,savefolderHRVtime, Neonate, win,Sessions(S,1).name,S)   
             disp('- Mean linelength finished')
        SDLL=SDLL_F(ECG_win_30,t_30,faktor,win);SavingC(SDLL,savefolderHRVtime, Neonate, win,Sessions(S,1).name,S) %Standart derivation of 5min linelength
             disp('- SDLL finsihed')
        SDaLL=SDaLL_F(ECG_win_30,t_30,faktor,win);SavingC(SDaLL,savefolderHRVtime, Neonate, win,Sessions(S,1).name,S) %Standart derivation of 30s linelength meaned over 5min
             disp('- SDaLL finished') 

  %%%%%%%% HRV TIME DOMAIN
        disp('HRV time domain analysis start')

        SDNN=SDNN_F(RR_300);SavingC(SDNN,savefolderHRVtime, Neonate,win,Sessions(S,1).name,S)
            disp('- SDNN finished') 
        RMSSD=RMSSD_F(RR_300);SavingC(RMSSD,savefolderHRVtime, Neonate,win,Sessions(S,1).name,S)
            disp('- RMSSD finished')  
        [NN50,NN30,NN20,NN10]=NNx(RR_300);SavingC(NN50,savefolderHRVtime, Neonate,win,Sessions(S,1).name,S);SavingC(NN30,savefolderHRVtime, Neonate,win,Sessions(S,1).name,S);SavingC(NN20,savefolderHRVtime, Neonate,win,Sessions(S,1).name,S);SavingC(NN10,savefolderHRVtime, Neonate,win,Sessions(S,1).name,S)
            disp('- NNx finished') 
        [pNN50,pNN30,pNN20,pNN10]=pNNx(RR_300);SavingC(pNN50,savefolderHRVtime, Neonate,win,Sessions(S,1).name,S);SavingC(pNN30,savefolderHRVtime, Neonate,win,Sessions(S,1).name,S);SavingC(pNN20,savefolderHRVtime, Neonate,win,Sessions(S,1).name,S);SavingC(pNN10,savefolderHRVtime, Neonate,win,Sessions(S,1).name,S)
            disp('- pNNx finished') 
        SDANN=SDANN_F(RR_30,faktor,win);SavingC(SDANN,savefolderHRVtime, Neonate,win,Sessions(S,1).name,S)
            disp('- SDANN finished')
        pDec=pDec_F(RR_300);SavingC(pDec,savefolderHRVtime, Neonate,win,Sessions(S,1).name,S)
            disp('- pDEC finished') 
        SDDec=SDDec_F(RR_300);SavingC(SDDec,savefolderHRVtime, Neonate,win,Sessions(S,1).name,S)
           disp('- SDDec finished')
           
  %%%%%%% HRV Frequency domain
        disp('Frequency time domain start')

        [totpow,VLF,LF,LFnorm,HF,HFnorm,ratioLFHF,sHF,sHFnorm,uHF,uHFnorm]=...
        freqdomainHRV (powerspectrum,f);SavingC(totpow,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S);SavingC(VLF,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S);SavingC(LF,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S);SavingC(LFnorm,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S);SavingC(HF,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S);SavingC(HFnorm,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S);SavingC(ratioLFHF,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S);SavingC(sHF,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S);SavingC(sHFnorm,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S);SavingC(uHF,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S);SavingC(uHFnorm,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S)
           disp('- Frequency finished') 
        [totpowR,LFR,LFnormR,HFR,HFnormR,ratioLFHFR,MFR,MFnormR,ratioMFHFR]=...
        freqdomainEDR (powerspectrumEDR,fEDR); SavingC(totpowR,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S); SavingC(LFR,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S);SavingC(LFnormR,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S); SavingC(HFR,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S); SavingC(HFnormR,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S); SavingC(ratioLFHFR,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S); SavingC(MFR,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S);   SavingC(MFnormR,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S);    SavingC(ratioMFHFR,savefolderHRVfreq, Neonate,win,Sessions(S,1).name,S);
           disp('- EDR requency finished')           

    %%%%%%% HRV Non linear
        disp('Nonlinear analysis start')
 
        [SampEn,QSE,SEAUC,r_opt]=SampEn_QSE_SEAUC(RR_300,faktor); SavingC(SampEn,savefolderHRVnonlin, Neonate, win,Sessions(S,1).name,S);SavingC(QSE,savefolderHRVnonlin, Neonate, win,Sessions(S,1).name,S);SavingC(SEAUC,savefolderHRVnonlin, Neonate, win,Sessions(S,1).name,S);SavingC(r_opt,savefolderHRVnonlin, Neonate, win,Sessions(S,1).name,S)
           disp('- SampEn QSE SEAUC finished')
        LZECG=LempelZivECG(ECG_win_300);  SavingC(LZECG,savefolderHRVnonlin, Neonate, win,Sessions(S,1).name,S)
          disp('- LepelZiv ECG finished')         
        LZNN=LempelZivRR(RR_300); SavingC(LZNN,savefolderHRVnonlin, Neonate, win,Sessions(S,1).name,S)
          disp('- LepelZiv HRV finished')   
        [SampEn_EDR,QSE_EDR,SEAUC_EDR,r_opt_EDR]=SampEn_QSE_SEAUC(EDR_300,faktor); SavingC(SampEn_EDR,savefolderHRVnonlin, Neonate, win,Sessions(S,1).name,S);SavingC(QSE_EDR,savefolderHRVnonlin, Neonate, win,Sessions(S,1).name,S);SavingC(SEAUC_EDR,savefolderHRVnonlin, Neonate, win,Sessions(S,1).name,S);SavingC(r_opt_EDR,savefolderHRVnonlin, Neonate, win,Sessions(S,1).name,S)
            disp('- SampEn_EDR QSE_EDR SEAUC_EDR finished')          
        LZEDR=LempelZivEDR(EDR_300); SavingC(LZEDR,savefolderHRVnonlin, Neonate, win,Sessions(S,1).name,S)
          disp('- LepelZiv EDR finished')  

        clearvars ECG_win_300 ECG_win_30 t_ECG_300 t_ECG_30 RR_idx_300 RR_300 RR_idx_30 RR_30 powerspectrum f powerspectrumEDR fEDR ECG Resp EMG EOG Chin...
            EDR_300 EDR_30 t_30 t_300        
     end %Sessionp
 end% Patient
 disp('----------------------------------')
 t1 = toc;
 dur=datestr(t1/(24*60*60),'DD:HH:MM:SS');
 disp('Finished' )
 disp (['Duration: ' dur]);  


 if shutdown
     pause('on');
     pause(20);
     system('shutdown -s')
 end
 
 %% Nested saving
    function Saving(Feature,savefolder, Neonate, win)
        if exist('Feature','var')==1
            name=inputname(1); % variable name of function input
            save([savefolder name '_win_' num2str(win) '_' num2str(Neonate)],'Feature')
        else
            disp(['saving of ' name ' not possible'])
        end       
    end
    function SavingAnnotations(Annotations,savefolder, Neonate, win)
        if exist('Annotations','var')==1
            name=inputname(1); % variable name of function input
            save([savefolder name '_win_' num2str(win) '_' num2str(Neonate)],'Annotations')
        else
            disp(['saving of ' name ' not possible'])
        end       
    end
    function SavingC(Feature,savefolder, Neonate, win,Session,S)
        if exist('Feature','var')==1
            name=inputname(1); % variable name of function input
            save([savefolder name '_Session_' num2str(S) '_win_' num2str(win) '_' Session],'Feature')
        else
            disp(['saving of ' name ' not possible'])
        end       
    end
 