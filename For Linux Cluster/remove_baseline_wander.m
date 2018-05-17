function [ECG_no_baselinewander]=remove_baseline_wander(ECG, Fs)

 %Remove basline wander with Butterworth bandpass 


%____________Create Filter_______________________________
Fn = Fs/2;                                          % Nyquist Frequency
Wp = [1  100]/Fn;                                   % Normalised Passband
Ws = [0.5  120]/Fn;                                 % Normalised Stopband
Rp = 10;                                            % Passband Ripple (dB)
Rs = 30;                                            % Stopband Ripple (dB)
[n,Ws] = cheb2ord(Wp, Ws, Rp, Rs);                  % Chebyshev Type II Order
[b,a] = cheby2(n, Rs, Ws);                          % Transfer Function Coefficients
[sos,g] = tf2sos(b,a);                              % Second-Order-Section For Stability
%_________________________________________________________
ECG_no_baselinewander = filtfilt(sos,g,ECG);

end
    
    
        
   
