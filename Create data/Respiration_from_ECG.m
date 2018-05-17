function [EDR]=Respiration_from_ECG(ECG,FS)
%----------------------------------------------------------------------
% ECG Derived Respiration (EDR)
% the resiration can be determined in seveeral ways:

% 1# Determine the R peaks-> Respiration signal->  interpolate to get same length as ECG signal
       % This is more a mechanical signal due to the movement of the chest and changing impedance
% 2# Determine HRV singal -> respiration signal 
        % ths is more a physiological signal which is the coupling between
        % many influences on the ECG via the respiration. 
% #3 Create a band pass filter -> remove everything exept low frequency parts -> REspiration signla
        % this is actually the same as the R peak detection. 
        
% the signal can be detrended beforehand to remove the very low frequency
% part. But be aware, maybe baselinedrift s also due to respiration...!

% AS IT IS KNOWN THAT THE RESPIRATION SINUS ARRYTHMIA IS NOT STRONG IN
% PRETERM INFANTS, WE RECCOMEND TO NLY USE METHOD 1 OR 3. MAYBE BASELINE
% DRIFT REMOVED. 

% The EDR of the intelleview signal is quite clear while not for the bin
% data. 

%Jan Werth

%----------------------------------------------------------------------
%% ******** ECG R peak detection **********
cd('C:\Users\310122653\Documents\PhD\cECG Data\Matlab\R peak detection')

% FS=500;
t_ecg=linspace(0,ceil(length(ECG)/FS), length(ECG))'; % timeline in seconds with 500 Hz FS    
hr_max=230;ploting=0;saving=0;

[ecg_r_peak_idx, ecg_s_peak_idx, ~, ~, ~, bbi_ecg, ~] = ecg_find_rpeaks(t_ecg, ECG, FS, hr_max,ploting,saving);

%% ******* 1# EDR from R peak amplitude not detrended ************

EDR_rpeak=ECG(ecg_r_peak_idx);                                                               % Respiration derived from ECG
EDR_FS=length(ecg_r_peak_idx)/(length(ECG)/FS);                                              % determine sampfle frequency for EDR signal to interpolate to Respiration signal (e.g.: bpm=120 -> sf=2Hz)
t_edr=linspace(0,floor(length(ecg_r_peak_idx)/EDR_FS), length(ecg_r_peak_idx))';             % Timeline in seconds with around 2 Hz fs    
EDR_same_length_as_ECG=interp1(t_edr,EDR_rpeak,t_ecg,'pchip'); 


%% ******* 2# EDR from HRV signal *******

EDR_HRV=bbi_ecg;

%% ******* 3# EDRform band pass filter *******

%____________Create Filter_______________________________
Fs = 500;                                           % Sampling Frequency
Fn = Fs/2;                                          % Nyquist Frequency
Wp = [0.2  2]/Fn;                                   % Normalised Passband
Ws = [0.1  3]/Fn;                                   % Normalised Stopband
Rp = 10;                                            % Passband Ripple (dB)
Rs = 30;                                            % Stopband Ripple (dB)
[n,Ws] = cheb2ord(Wp, Ws, Rp, Rs);                  % Chebyshev Type II Order
[b,a] = cheby2(n, Rs, Ws);                          % Transfer Function Coefficients
[sos,g] = tf2sos(b,a);                              % Second-Order-Section For Stability
%_________________________________________________________
EDR_passband = filtfilt(sos,g,ECG); %remove baseline wander with zero phase filter

%% *********************
%******* Decide which variable to get out *******
EDR=EDR_same_length_as_ECG;
EDR=EDR-nanmean(EDR); %center around zero. Remove offset.

%%

% %Remove outliers in ECG due to lead switch
% EDR_passband(EDR_passband>mean(EDR_passband)+6*std(EDR_passband)) =[]; EDR_passband(EDR_passband<mean(EDR_passband)-6*std(EDR_passband))=[];
% EDR_resp_500(EDR_resp_500>mean(EDR_resp_500)+6*std(EDR_resp_500)) =[]; EDR_resp_500(EDR_resp_500<mean(EDR_resp_500)-6*std(EDR_resp_500))=[];
% 
% 
%  %Getting Respiration data on same level
%  Resp_intelleview_data_tot_500_scaled = -1 + 2.*(Resp_intelleview_data_tot_500 - min(Resp_intelleview_data_tot_500))./(max(Resp_intelleview_data_tot_500) - min(Resp_intelleview_data_tot_500));
%  EDR_passband_scaled = -1 + 2.*(EDR_passband - min(EDR_passband))./(max(EDR_passband) - min(EDR_passband));
%  EDR_resp_500_scaled = -1 + 2.*(EDR_resp_500 - min(EDR_resp_500))./(max(EDR_resp_500) - min(EDR_resp_500));



end