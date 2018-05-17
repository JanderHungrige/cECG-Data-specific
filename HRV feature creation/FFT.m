%Frequency analysis
% for un-equidistand data frequency analysis look at:
% http://nl.mathworks.com/matlabcentral/newsreader/view_thread/41005

% first the data has to be interpoliated to be equidistand in time
% then the FFT can be calculated to generate the spectrum

% clc
% clear
% win=[60,120,180,240,300]; %window in seconds
% overlap=[0.5,1];

function FFT(pat,plotting)
    pat=[3];
    for i=1:length(pat)
                    cd('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\ECG no care separated into AS and QS (total)\HRV analysis\Ralphs\')
%%%%%%%%%%%%%%%%%AS                    
                   if exist(fullfile(cd, ['RRdistanceAS_total_' num2str(pat(1,i)) '.mat']), 'file') == 2 % ==> 0 or 2
                        load(fullfile(cd, ['RRdistanceAS_total_' num2str(pat(1,i)) '.mat']))                       
                        load(['C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\ECG no care separated into AS and QS (total)\ECGstate_nocare_' num2str(pat)])%ECG file for timevector creation
                        RRdistanceAS(any(isnan(RRdistanceAS)))=[]; %removing nans
                       
                       

                        timing=(AS(1,1):length(RRdistanceAS):AS(1,end)); % create time vector from ECG file
%                         timing=cumsum(RRdistanceAS);%creating increasing time vector used for t
                        
                        %windowing
                        
                        %interpolation
                        
                        
                        %FFT
                        
                        NFFT = 2^nextpow2(RRdistanceAS);    
                        Y = fft(y,NFFT)/RRdistanceAS;
                        f = linspace(0, Fs/2, NFFT/(2+1));
                        
                        %Using FFT to learn more about frequency components present in your signal is RRdistance. 
                        %You can use the function pwelch function in MATLAB to learn more frequencies present in your 
                        %signal and also the power of these signals. MATLAB will automatically compute the NFFT required 
                        %and return the frequencies present in your signal along with the power at each frequency. Use this syntax: [p,f] = pwelch(x,[],[],[],Fs)

                        %Look at the documentation of pwelch for more information.
                        
                        
                        
                        %plotting
                        if plotting
                            figure
                            plot(f, Y)
                            title (['Powerspectrum AS patient: ' num2str(pat)])
                        end
                        
                        
                   end
%%%%%%%%%%%%%%%%QS
                    if exist(fullfile(cd, ['RRdistanceQS_total_' num2str(pat(1,i)) '.mat']), 'file') == 2 % ==> 0 or 2
                        load(fullfile(cd, ['RRdistanceQS_total_' num2str(pat(1,i)) '.mat']))                       
                        
                        RRdistanceQS(any(isnan(RRdistanceQS)))=[]; %removing nans


                        timing=cumsum(RRdistanceQS);%creating increasing time vector used for t
                        
                        %windowing
                        
                        %interpolation
                        
                        
                        %FFT
                        
                        NFFT = 2^nextpow2(RRdistanceQS);    
                        Y2 = fft(y,NFFT)/RRdistanceQS;
                        f2 = linspace(0, Fs/2, NFFT/2+1));
                        
                        %plotting
                        if plotting
                            figure
                            plot(f2, Y2)
                            title (['Powerspectrum QS patient: ' num2str(pat)])
                        end
                                               
                    end
    end
    
    
    
end
