    %afterwards
function [signal_filtered]=hilbertdemod(x,filterCoefficients,first,overlap)%[signal_filtered]=hilbertdemod(x,fc_lpf,fs,fc_lpf2,fc_hpf, filterCoefficients,first)
%inputs are: x: signal, fc_lpf: low pass filter cut off (1 Khz), fs:
%sampling frequency. fc_lpf2; second low pass filter stage.
%outputs are: signal_filtered, the filtered signal.
%an example is:     signal_filtered=hilbertdemod(signal(1:maxpt),1000,8000,100);
%s=size(x);
%sample_size=s(:,2);
%t=(1/fs:1/fs:sample_size/fs);
%High pass filter the signal

persistent prev;
persistent prev1;
persistent prev2;

%High pass 100
if isempty(first)
    [prev,x1]=iirFilter([],x',filterCoefficients.iir_high,'first',overlap);
else
    [prev,x1]=iirFilter(prev,x',filterCoefficients.iir_high,'',overlap);
end

%    fNorm = 100 / (fs/2);
%    [d,c] = butter(8, fNorm, 'high');
%    x= filtfilt(d, c, x);
   
   
signal=abs(hilbert(x1));

% rsig1=sin(2*pi*1000*t);
% rsig2=cos(2*pi*1000*t);

%Low pass 1000
if isempty(first)
    [prev1,x2]=iirFilter([],signal,filterCoefficients.iir_low1,'first',overlap);
else
    [prev1,x2]=iirFilter(prev1,signal,filterCoefficients.iir_low1,'',overlap);
end
% fNorm = 1000 / (fs/2);
% [b,a] = butter(8, fNorm, 'low');
% signal_filtered= filtfilt(b, a,(signal1));


%Low pass 100
if isempty(first)
    [prev2,signal_filtered]=iirFilter([],x2,filterCoefficients.iir_low2,'first',overlap);
else
    [prev2,signal_filtered]=iirFilter(prev2,x2,filterCoefficients.iir_low2,'',overlap);
end
% fNorm = 100 / (fs/2);
% [b1,a1] = butter(8, fNorm, 'low');
% signal_filtered= filtfilt(b1,a1,(signal_filtered));


end
   
