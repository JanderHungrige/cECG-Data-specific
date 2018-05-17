%1)This function loads rspiration data and finds via peak detetion and peak
%width if there is clipping.
% the peak detecktion is limited to plateau peaks -> clipped peaks.
%2) Later the inteleview ECG is loaded and the respiration is derived from
% there. 
%3)The ECG and Resp are synchronized
%4) The clipped resiration is completed with the ECG derived Respiration

clear
clc


%-----------------------------------
%Load Respiration data
%-----------------------------------
dataSetpath = ('\\code1\storage\2012-0194_neonatal_data\cECG study\20120719_test9\Intelliview serial datalogging\2012-07-19 12-33_001'); 
cd(dataSetpath)
Resp_intelleview_files=dir('Resp_*.txt'); % find all REf ECG files
numFiles=size(Resp_intelleview_files,1);%amount of ECG ref files 
for fileNumber=1:numFiles

filename=Resp_intelleview_files(fileNumber).name;
delimiter = '\t';
startRow = 15;

formatSpec = '%*s%f%[^\n\r]';
fileID = fopen(filename,'r');
textscan(fileID, '%[^\n\r]', startRow-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'ReturnOnError', false);
fclose(fileID);
Resp_intelleview_data{fileNumber} = [dataArray{1:end-1}];
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

end
cellsz = cellfun(@sum,cellfun(@numel,Resp_intelleview_data,'uni',false),'uni',false);
Resp_inteleview_length=sum(cellfun(@sum,cellsz)); clearvars cellsz;
    
       %reate one long Respiration file
Resp_intelleview_data_tot = Resp_intelleview_data{1,1};
if numFiles > 1
    for j=2:numFiles
    Resp_intelleview_data_tot=[Resp_intelleview_data_tot; Resp_intelleview_data{1,j}];
    end
end
    
clearvars j numFiles fileNumber
    
%-----------------------------------
%Load ECG data
%-----------------------------------
    
ECG_intelleview_files=dir('II_*.txt'); % find all REf ECG files
kb=ECG_intelleview_files.bytes;
if kb <= 1000 %data smaller than 1kbyte goto II...txt
    ECG_intelleview_files=dir('I_*.txt'); kb=ECG_intelleview_files.bytes;
    if kb <= 1000 %data still smaller than 1kbyte goto III...txt
        ECG_intelleview_files=dir('III_*.txt'); 
    end
end; clearvars kb
        
numFiles=size(ECG_intelleview_files,1);%amount of ECG ref files 
for fileNumber=1:numFiles
%read in txt files
filename=ECG_intelleview_files(fileNumber).name;
delimiter = '\t';
startRow = 15; % first few lines are nan. Check if always 15

formatSpec = '%*s%f%[^\n\r]';
fileID = fopen(filename,'r');
textscan(fileID, '%[^\n\r]', startRow-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'ReturnOnError', false);
fclose(fileID);
ECG_intelleview_data_cell{fileNumber} = [dataArray{1:end-1}];
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
end
% calc total length of intelleview ECG
cellsz = cellfun(@sum,cellfun(@numel,ECG_intelleview_data_cell,'uni',false),'uni',false);
ECG_inteleview_length=sum(cellfun(@sum,cellsz)); clearvars cellsz;

   %reate one long ECG file
ECG_intelleview_data_tot = ECG_intelleview_data_cell{1,1};
if numFiles > 1
    for j=2:numFiles
    ECG_intelleview_data_tot=[ECG_intelleview_data_tot; ECG_intelleview_data_cell{1,j}];
    end
end; clearvars numFiles j fileNumber


%Remove basline wander with Butterworth bandpass 
%____________Create Filter_______________________________
Fs = 500;                                           % Sampling Frequency
Fn = Fs/2;                                          % Nyquist Frequency
Wp = [1  100]/Fn;                                   % Normalised Passband
Ws = [0.5  120]/Fn;                                 % Normalised Stopband
Rp = 10;                                            % Passband Ripple (dB)
Rs = 30;                                            % Stopband Ripple (dB)
[n,Ws] = cheb2ord(Wp, Ws, Rp, Rs);                  % Chebyshev Type II Order
[b,a] = cheby2(n, Rs, Ws);                          % Transfer Function Coefficients
[sos,g] = tf2sos(b,a);                              % Second-Order-Section For Stability
%_________________________________________________________
ECG_detrend = filtfilt(sos,g,ECG_intelleview_data_tot); %remove baseline wander 

%-----------------------------------
%Sync Respiration and ECG
%-----------------------------------
Resp_FS=Resp_inteleview_length/(ECG_inteleview_length/500);
FS_factor=500/Resp_FS;
t_62=linspace(0,ceil(Resp_inteleview_length/Resp_FS), Resp_inteleview_length)'; % timeline in seconds with 62.5 Hz FS
t_500=linspace(0,ceil(Resp_inteleview_length/Resp_FS), Resp_inteleview_length*FS_factor)'; % timeline in seconds with 500Hz FS
Resp_intelleview_data_tot_500=interp1(t_62,Resp_intelleview_data_tot,t_500,'pchip');%Shape-preserving piecewise cubic interpolation

clearvars t_62 t_500
 

%----------------------------------------------------------------------------
% Repsiration Peak detection and plateau determination for Respiration signal
%----------------------------------------------------------------------------
maxmax=max(Resp_intelleview_data_tot_500);minmax=max(-Resp_intelleview_data_tot_500);
[pks_max,locs_max,w_max,p_max]=findpeaks(Resp_intelleview_data_tot,'MinPeakHeight',maxmax-0.001);
[pks_min,locs_min,w_min,p_min]=findpeaks(-Resp_intelleview_data_tot,'MinPeakHeight',minmax-0.001);

% check if there is clipping by mean length of peaks greater than 20
% samples. MEaning there is are peak plateaus. 
if  mean(w_max) >=20
    clipped=1;
else
    clipped=0;
    disp('THe Respiration data seems not to be clipped. Please check manually')
end

clearvars   maxmax   minmax pks_max pks_min p_max p_min

%----------------------------------------------------------------------
% ECG Derived Respiration (EDR)
%----------------------------------------------------------------------
%ECG R peak detection
fs=500; t_ecg=linspace(0,ceil(ECG_inteleview_length/fs), ECG_inteleview_length)'; % timeline in seconds with 500 Hz FS    
hr_max=230;ploting=0;saving=0;
cd('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\R peak detection and HRV\Ralps Rpeak detection')
[ecg_r_peak_idx, ~, ~, ~, ~, ~, ~] = ecg_find_rpeaks(t_ecg, ECG_detrend, fs, hr_max,ploting,saving);

%EDR determination via R peak amplitude
EDR=ECG_detrend(ecg_r_peak_idx);                                                             % Respiration derived from ECG
EDR_FS=length(ecg_r_peak_idx)/(ECG_inteleview_length/500);                                   % determine sampfle frequency for EDR signal to interpolate to Respiration signal (e.g.: bpm=120 -> sf=2Hz)
Resp_FS=Resp_inteleview_length/(ECG_inteleview_length/500);                                  % Determine sample Frequency for Respiration signal (should be 62.5 Hz)
t_edr=linspace(0,floor(length(ecg_r_peak_idx)/EDR_FS), length(ecg_r_peak_idx))';             % Timeline in seconds with around 2 Hz fs    
FS_factor=Resp_FS/EDR_FS;                                                                    % Determine the upsamppling factor from around 2 Hz form EDR to 62.5 Hzfrom Respiration
t_resp=linspace(0,ceil(length(ecg_r_peak_idx)/EDR_FS), length(ecg_r_peak_idx)*FS_factor)';   % Timeline in seconds with 60.2 Hz FS
EDR_resp=interp1(t_edr,EDR,t_resp,'pchip');     
EDR_resp_500=interp1(t_edr,EDR,t_ecg,'pchip');                                                  % Shape-preserving piecewise cubic interpolation of EDR signal to match length of Respiration 

%EDR determination via frequency band band pass filtering
%____________Create Filter_______________________________
Fs = 500;                                           % Sampling Frequency
Fn = Fs/2;                                          % Nyquist Frequency
Wp = [0.2  2]/Fn;                                 % Normalised Passband
Ws = [0.1  3]/Fn;                                   % Normalised Stopband
Rp = 10;                                            % Passband Ripple (dB)
Rs = 30;                                            % Stopband Ripple (dB)
[n,Ws] = cheb2ord(Wp, Ws, Rp, Rs);                  % Chebyshev Type II Order
[b,a] = cheby2(n, Rs, Ws);                          % Transfer Function Coefficients
[sos,g] = tf2sos(b,a);                              % Second-Order-Section For Stability
%_________________________________________________________
EDR_passband = filtfilt(sos,g,ECG_intelleview_data_tot); %remove baseline wander 



%Getting Respiration data on same level
%zero mean
Resp_intelleview_data_tot_500=Resp_intelleview_data_tot_500-mean(Resp_intelleview_data_tot_500(:));
EDR_passband=EDR_passband-mean(EDR_passband(:));
EDR_resp_500=EDR_resp_500-mean(EDR_resp_500(:));

figure
plot(Resp_intelleview_data_tot_500)
hold on
plot(EDR_passband*5,'r')
plot(EDR_resp_500*5,'c')
legend('Impedance','EDR-RR peaks','EDR-band pass')

 %Remove outliers in ECG due to lead switch
 EDR_passband(EDR_passband>mean(EDR_passband)+6*std(EDR_passband)) =[]; EDR_passband(EDR_passband<mean(EDR_passband)-6*std(EDR_passband))=[];
 EDR_resp_500(EDR_resp_500>mean(EDR_resp_500)+6*std(EDR_resp_500)) =[]; EDR_resp_500(EDR_resp_500<mean(EDR_resp_500)-6*std(EDR_resp_500))=[];


 %Getting Respiration data on same level
 Resp_intelleview_data_tot_500_scaled = -1 + 2.*(Resp_intelleview_data_tot_500 - min(Resp_intelleview_data_tot_500))./(max(Resp_intelleview_data_tot_500) - min(Resp_intelleview_data_tot_500));
 EDR_passband_scaled = -1 + 2.*(EDR_passband - min(EDR_passband))./(max(EDR_passband) - min(EDR_passband));
 EDR_resp_500_scaled = -1 + 2.*(EDR_resp_500 - min(EDR_resp_500))./(max(EDR_resp_500) - min(EDR_resp_500));
 
figure
plot(Resp_intelleview_data_tot_500_scaled)
hold on
plot(EDR_passband_scaled,'r')
plot(EDR_resp_500_scaled,'c')
legend('Impedance','EDR-RR peaks','EDR-band pass')

figure
hold on
plot(t_resp,EDR_resp,'r')
plot(t_62,Resp_intelleview_data_tot)
hold off
figure
hold on
plot(Resp_intelleview_data_tot)
plot(locs_max,Resp_intelleview_data_tot(locs_max),'rv','MarkerFaceColor','r')
plot(locs_min,Resp_intelleview_data_tot(locs_min),'rs','MarkerFaceColor','b')    


