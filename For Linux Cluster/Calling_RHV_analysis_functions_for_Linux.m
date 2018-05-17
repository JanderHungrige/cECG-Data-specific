function Calling_RHV_analysis_functions_for_Linux (varargin)
clear
clc
tic


patient=varargin{1};
loadfolder=varargin{2};
savefolder=varargin{3};
savefolderAnnotations=varargin{4};
savefolderHRVtime=varargin{5};
savefolderHRVfreq=varargin{6};
savefolderHRVnonlin=varargin{7};

pat=str2double(patient);

saving=1;
plotting=0;
win=300;
winXi=300;
faktor=30; % how much is the data moving forward? 30s is classic
FS_ecg=500;loadfolder

 for I=1:length(pat)
    disp('***************************************')
    disp(['Working on patient ' num2str(pat(I))])
    Neonate=pat(I);      
    
  
%% ************ Load data **************

    filetype='Intellivue';
    filelocation=([loadfolder 'participant' num2str(Neonate) '/' filetype '/']);
    Datum=dir([filelocation '2012*']); % folder name differnt on this stage
    filelocation=([filelocation Datum.name '/']);
    Sessions=dir([filelocation filetype '*']);
    
    for S=1:length(Sessions)
        disp('- - - - - - - - - - - - - - - - - - - ')
        disp(['Working on session: ' Sessions(S,1).name])
        load([filelocation Sessions(S,1).name])
        
        
%% ************ Load annotations (1s) **************  
%loading 1 secondannotations for this particular patient/session

    Annotation=loading_annotations(Neonate,Sessions(S,1).name);
    disp('* Annotation loaded')
                
 %% ************ Window  ECG /  Annotation signals 


    t_ECG=linspace(0,length(ECG.values)/FS_ecg,length(ECG.values))';
    % The differnec in t_300 and t_ECG_300 is that t_ECG_100 is a
    % continuous run of time, while t_300 is 0 to t for each cell element

   [ECG_win_300,ECG_win_30,t_ECG_300,t_ECG_30]=SlidingWindow_ECG(ECG.values,t_ECG,Neonate,win,saving,savefolder,faktor); 
   [Annotations_win_300, Annotations_win_30]=SlidingWindow_Annotations(Annotation,t_ECG,Neonate,win,saving,savefolderAnnotations,faktor);
     disp(['* Data is merged into windows of length: ' num2str(win) 's and ' num2str(30) 's'] )        
%% ************ Creating RR signal for ECG-Signal **************
    Ralphsfactor={1;-1;-1;1;1;-1;1;-1; 1; 1;-1; 1;-1; 1;-1;-1;-1;-1};%Determine if the ECG signal should be turned -1 or not 1. 
                 %1  2 3  4 5  6 7  8  9  10 11 12 13 14 15 16 17 18
    padding=0; %Determine if the RR should be same length as ECG. Don`t have to be
    plotting=0; %plotting Ralphs RR detection

       
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
         
%% ************ Creating spectrum for ECG-Signal **************         
   [powerspectrum,f]=Lomb_scargel_single(RR_300,RR_idx_300,t_300,Neonate,saving,savefolderHRVtime,win,Sessions(S,1).name,S) ;
    disp('* Periodogram calculated')
         


%% ************ calculate HRV **************

%%%%%%%% ECG TIME DOMAIN     
      disp('ECG time domain analysis start')        

    Beats_per_Epoch(RR_300,Neonate,saving,savefolderHRVtime,win,Sessions(S,1).name,S)   % S for session number
         disp('- BpE finished')
        
    linelength(ECG_win_300,t_300,Neonate,saving,savefolderHRVtime,win,Sessions(S,1).name,S)     
         disp('- Linelength finished')
    meanarclength(ECG_win_30,t_30,Neonate,saving,savefolderHRVtime,win,faktor,Sessions(S,1).name,S) 
         disp('- Mean linelength finished')
    SDLL(ECG_win_30,t_30,Neonate,saving,savefolderHRVtime,win,faktor,Sessions(S,1).name,S) %Standart derivation of 5min linelength
         disp('- SDLL finsihed')
    SDaLL(ECG_win_30,t_30,Neonate,saving,savefolderHRVtime,win,faktor,Sessions(S,1).name,S) %Standart derivation of 30s linelength meaned over 5min
         disp('- SDaLL finished')


%%%%%%%% HRV TIME DOMAIN
    disp('HRV time domain analysis start')
 
    SDNN(RR_300,Neonate,saving,savefolderHRVtime,win,Sessions(S,1).name,S);
        disp('- SDNN finished') 
    RMSSD(RR_300,Neonate,saving,savefolderHRVtime,win,Sessions(S,1).name,S);
        disp('- RMSSD finished')  
    NNx(RR_300,Neonate,saving,savefolderHRVtime,win,Sessions(S,1).name,S);
        disp('- NNx finished') 
    pNNx(RR_300,Neonate,saving,savefolderHRVtime,win,Sessions(S,1).name,S);
        disp('- pNNx finished') 
    SDANN(RR_30,Neonate,saving,savefolderHRVtime,win,faktor,Sessions(S,1).name,S);
        disp('- SDANN finished')  

    
%%%%%%% HRV Frequency domain

    disp('Frequency time domain start')
    freqdomainHRV (powerspectrum,f,Neonate,win,saving,savefolderHRVfreq,Sessions(S,1).name,S)

%%%%%%% HRV Non linear
    disp('nonlinearr analysis start')
    
    SampEn_QSE(RR_300,Neonate,saving,savefolderHRVnonlin,win,faktor,Sessions(S,1).name,S ) %
    disp('- SampEn and QSE finished')


    clearvars ECG_win_300 ECG_win_30 t_ECG_300 t_ECG_30 RR_idx_300 RR_300 RR_idx_30 RR_30 powerspectrum f
    end %Session

 end% Patient
toc
end