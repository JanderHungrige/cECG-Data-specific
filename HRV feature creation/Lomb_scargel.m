%FRequency analysis

%[pxx,f] = plomb(x,t) returns the Lomb-Scargle power spectral density (PSD) estimate, pxx, of a signal, x, that is sampled at the instants specified in t. 
% t must increase monotonically but need not be uniformly spaced. All elements of t must be nonnegative. pxx is evaluated at the frequencies returned in f.

clc
clear
win=[60,120,180,240,300]; %window in seconds
overlap=[0.5,1];
function Lomb_scargel(win, overlap,pat,plotting)
    pat=[3];
    for i=1:length(pat)
        for j=1%:length(win)
            for o=1%:length(overlap)
                    cd('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\windowed_ECG_chunks\Ralphs\')
                   if exist(fullfile(cd, ['RRdistanceAS' num2str(pat(1,i)) ' _ ' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,o)) '.mat']), 'file') == 2 % ==> 0 or 2
                        load(fullfile(cd, ['RRdistanceAS' num2str(pat(1,i)) ' _ ' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,o)) '.mat']))                       
                        
                        timing{1,i}=cumsum(RRdistanceAS{1,i});%creating increasing time vector used for t
                        
                        %windowing
                        
                        
                        %lomb scargel
                        
                        %normalizing
                        
                        %adding single windows
                       
                   end
                   if exist(fullfile(cd, ['RRdistanceQS' num2str(pat(1,i)) ' _ ' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,o)) '.mat']), 'file') == 2 % ==> 0 or 2

                         load(fullfile(cd, ['RRdistanceQS' num2str(pat(1,i)) ' _ ' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,o)) '.mat']))
                   
                        timing{1,i}=cumsum(RRdistanceQS{1,i});%creating time vector used for t
                         
                        %windowing
                        
                        
                        %lomb scargel
                        
                        %normalizing
                        
                        %adding single windows
                   end






            end
        end
    end
end
