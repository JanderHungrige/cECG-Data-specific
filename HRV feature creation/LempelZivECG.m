function LempelZivECG(ECG,Neonate,saving,savefolder,win,Session,S)  
%Input
% RR: 5min RR distance data
% Neonate: Which patient
% saving: If saving is whished
% savefolder: Where to save
% win: Duration of the HRV window. Comon is 5min/300s

% One way to do it would be to take the Hilbert transform (matlab function 'hilbert'),
% of the timeseries data to get the analytic signal. 
% Taking the absolute value of this (matlab function 'abs'), 
% you get a new set of time series indicating the amplitude of the signal
% at each timepoint. You can binarize the data by thresholding this this time series.
% Schartner et al 2015-2016 suggest using the mean of the signal amplitude as a
% threshold. ref: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0133532

% https://nl.mathworks.com/matlabcentral/fileexchange/38211-calc-lz-complexity

%   The type of complexity to evaluate as a string, which is one of:
%       - 'exhaustive': complexity measurement is based on decomposing S 
%       into an exhaustive production process.
%       - 'primitive': complexity measurement is based on decomposing S 
%       into a primitive production process.
%   Exhaustive complexity can be considered a lower limit of the complexity
%   measurement approach proposed in LZ76, and primitive complexity an
%   upper limit.
%
%   normalize:
%   A logical value (true or false), used to specify whether or not the 
%   complexity value returned is normalized or not.  
%   Where normalization is applied, the normalized complexity is 
%   calculated from the un-normalized complexity, C_raw, as:
%       C = C_raw / (n / log2(n))
%   where n is the length of the sequence S.

type='exhaustive'; % exhaustive or primitive
normalize=1 ;% 1 or 0

% down sampling to 20Hz
for i=1:length(ECG)
    ECG{1,i}=downsample(ECG{1,i},500/20);
end


for i=1:length(ECG)
    if all(isnan(ECG{1,i}))
        LZECG{i,1}=nan;
        continue
    end
    ECG{1,i}(isnan(ECG{1,i}(:,1)))=[]; % this removes nans for the hilbert
    Hilberttimeseries=abs(hilbert(ECG{1,i}(:,1))); % we use 2:end to avoide the nan at the beginning
    thres=nanmean(Hilberttimeseries); % determine threshold forom the hilbert time series
    Binarieinz=find(Hilberttimeseries>thres);
    Binarinull=find(Hilberttimeseries<thres);
    Hilberttimeseries(Binarieinz,1)=1; 
    Hilberttimeseries(Binarinull,1)=0;
    Hilbertsrting = binary_seq_to_string(Hilberttimeseries);
    [LZECG{i,1},~,~] = calc_lz_complexity(Hilbertsrting, type, normalize);
end
    

%%%%%%%%%%%%replace [] with nan
ix=cellfun(@isempty,LZECG);
LZECG(ix)={nan};  

            
%%%%%%%%%%%% SAVING            
if saving                     %saving R peaks positions in mat file                 
    Saving(LZECG,savefolder,Neonate,win,Session,S) 
end% end if saving 

    
  
end

%% Nested saving
    function Saving(Feature,savefolder, Neonate, win,Session,S)
        if exist('Feature','var')==1
            name=inputname(1); % variable name of function input
            save([savefolder name '_Session_' num2str(S) '_win_' num2str(win) '_' Session],'Feature')
        else
            disp(['saving of ' name ' not possible'])
        end       
    end
 
 
 

