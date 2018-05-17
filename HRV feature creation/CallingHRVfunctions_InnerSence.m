%Calling RHV analysis functions
clear
clc
tic

Matlabbase='C:\Users\310122653\Documents\PhD\Matlab\InnerSence Data specific\HRV feature creation\';
addpath('C:\Users\310122653\Documents\PhD\Matlab');

path='E'; % the HDD file with the patient data. Needed for loading ECG

pat=[3,4,5,6,7,9,13,15]; 
% pat=5
saving=0;
plotting=0;
win=[30 60 150 300];
winXi=300;


 for i=1:length(pat)
    disp('***************************************')
    disp(['Working on patient ' num2str(pat(i))])
    Neonate=pat(i);        
        
    
% %%%%%%%% ECG TIME DOMAIN   
%         cd(Matlabbase)   
%         loadfolder=([path ':\InnerSence\processed data\ECG\']);
%         savefolderHRVtime= ([path ':\InnerSence\processed data\HRV-features\timedomain\']);
%     Beats_per_Epoch(Neonate,saving,win,loadfolder,savefolderHRVtime)
%     
% % cd('C:\Users\310122653\Documents\PhD\Matlab\ECG analysis')
% %         if exist('arclenAS','var')==0
% %             arclenAS=[];arclenQS=[];
% %         end
% %         cd('C:\Users\310122653\Documents\PhD\EEG_MMC Data\Matlab\ECG analysis')
% %         [arclenAS,arclenQS]=linelength(Neonate,saving,arclenAS,arclenQS);
% %         cd('C:\Users\310122653\Documents\PhD\EEG_MMC Data\Matlab\ECG analysis')
% %          meanarclength(Neonate,pat,arclenAS,arclenQS);
%
% 
% %%%%%%%% HRV TIME DOMAIN
%         disp('HRV time domain start')
%         loadfolder=([path ':\InnerSence\processed data\ECG\']);
%         savefolderHRVtime= ([path ':\InnerSence\processed data\HRV-features\timedomain\']);
%         cd(Matlabbase)   
%     SDNN(Neonate,saving,win,loadfolder,savefolderHRVtime);
% %         disp('SDNN finished') 
% %         cd(Matlabbase)
%     RMSSD(Neonate,saving,win,loadfolder,savefolderHRVtime);
% %         disp('RMSSD finished')  
% %         cd(Matlabbase)
%     NNx(Neonate,saving,win,loadfolder,savefolderHRVtime);
% %         disp('NNx finished') 
% %         cd(Matlabbase)
%     pNNx(Neonate,saving,win,loadfolder,savefolderHRVtime);
%         disp('pNNx finished') 
%         cd(Matlabbase) 
% %     SDANN(Neonate,saving,win,loadfolder,savefolderHRVtime);
% %         disp('SDANN finished')  
% %         cd(Matlabbase)        
%         disp('HRV time domain finished')   
% 
%     
% %%%%%%% HRV Frequency domain
%         disp('Frequency time domain start')
%         loadfolder=([path ':\InnerSence\processed data\ECG\Spectrum\']);
%         savefolderHRVfreq= ([path ':\InnerSence\processed data\HRV-features\freqdomain\']);        
%         cd(Matlabbase)
%    freqdomainHRV (Neonate,plotting,win,saving,loadfolder,savefolderHRVfreq)
%         cd(Matlabbase)
%  disp('Frequency time domain finished')   
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - XI - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(Matlabbase)  
Xifolder=genpath('C:\Users\310122653\Documents\PhD\Matlab\Xi features\ecg');
addpath(Xifolder);

settings.HRV_SPECTRAL_BANDS= [   0.04,0.15;...
                                 0.15,0.4;...
                                 0.4, 0.7 ];  %+[0.7, 1.5] three bands expected
settings.EPOCH_LENGTH=30; % in s
settings.HRV_WINDOW_LENGTH=300; % 5min
settings.CONSIDER_COVERAGE=1; % yes/no
settings.PROCESSING_THR_WINDOW_COVERAGE= 1;% yes/no

            % settings.HRV_WINDOW_LENGTH being expressed in seconds, assuming an
            % average of 1 beat per second, the computation of features is
            % performed only if enough RR-intervals are present

%%%%%%% loading RR data
        loadfolder=([path ':\InnerSence\processed data\ECG\']);

        [consence_RR_windowed,winlength]=load_RR(loadfolder,Neonate,winXi);
        if winlength==1
            return
        end; clearvars winlength

%%%%%%% Get going 

        loadfolder=([path ':\InnerSence\processed data\ECG\']);
        savefolderHRVtime= ([path ':\InnerSence\processed data\HRV-features\Xi\']);
        
ret = extract_ecg_hrv_fd_adpt(consence_RR_windowed, settings, savefolderHRVtime);


 end
toc