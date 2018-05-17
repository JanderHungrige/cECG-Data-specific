%Frequency analysis

%[pxx,f] = plomb(x,t) returns the Lomb-Scargle power spectral density (PSD) estimate, pxx, of a signal, x, that is sampled at the instants specified in t. 
% t must increase monotonically but need not be uniformly spaced. All elements of t must be nonnegative. pxx is evaluated at the frequencies returned in f.

% clc
% clear
% win=[60,120,180,240,300]; %window in seconds
% overlap=[0.5,1];
function Lomb_scargel_windowed(pat,plotting,win,overlap,saving)
    for i=1:length(pat)
        for j=1:length(win)
            for k=1:length(overlap)
                %Loading RR distance data form Ralphs RR distance calculations
                    cd('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\windowed_ECG_chunks\Ralphs\')
                   if exist(fullfile(cd, ['RRdistanceAS' num2str(pat(1,i)) ' _ ' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,k)) '.mat']), 'file') == 2 % ==> 0 or 2
                        load(fullfile(cd, ['RRdistanceAS' num2str(pat(1,i)) ' _ ' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,k)) '.mat'])) % load  RR differnces 
                        
                    %removing nans. Due to the RR distance calculation the first value in NAN
                        for l=1:length(RRdistanceAS)
                            RRdistanceAS{1,l}(any(isnan(RRdistanceAS{1,l})))=[]; %removing nans
                        end

                    % removing outlier
                        for l=1:length(RRdistanceAS)
                            if isempty(RRdistanceAS{1,l})==0
                                all_idx = 1:length(RRdistanceAS{1,l});
                                outlier_idx = abs(RRdistanceAS{1,l} - median(RRdistanceAS{1,l})) > 3*std(RRdistanceAS{1,l}); % Find outlier idx
                                RRdistanceAS{1,l}(outlier_idx) = interp1(all_idx(~outlier_idx), RRdistanceAS{1,l}(~outlier_idx), all_idx(outlier_idx)); % Linearly interpolate over outlier idx for x
                            end
                        end
                        
                 %Loading R peaks from Raphs R peak detection      
                       if exist(fullfile(cd, ['Rpeaks' num2str(pat(1,i)) ' _ ' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,k)) '.mat']), 'file') == 2 % ==> 0 or 2
                            load(fullfile(cd, ['Rpeaks' num2str(pat(1,i)) ' _ ' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,k)) '.mat'])) % load R peaks
                       end
                       
                 %Loading ECG and time for cration of timevector regarding the R peaks from Ralph      
                    cd('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\ECG_chunks')
                       if exist(fullfile(cd, ['ChunkAS' num2str(pat(1,i)) '_win' num2str(win(1,j)) '_overlap' num2str(overlap(1,k)) '.mat']), 'file') == 2 % ==> 0 or 2
                            load(fullfile(cd, ['ChunkAS' num2str(pat(1,i)) '_win' num2str(win(1,j)) '_overlap' num2str(overlap(1,k)) '.mat']))  % load QS and AS; data and time vector of all AS/QS
       
                           % creating time vectorfor AS
                           for l=1:length(RpeakAS)
                               if isempty(ChunkAS{1,l}) == 0 
                                  timingAS{1,l}=ChunkAS{1,l}(1,RpeakAS{1,l}(2:end)); % get the timestamp from AS time vector at the Rpeaks in s
                               end
                           end
                       end
                           
                       if exist(fullfile(cd, ['ChunkQS' num2str(pat(1,i)) '_win' num2str(win(1,j)) '_overlap' num2str(overlap(1,k)) '.mat']), 'file') == 2 % ==> 0 or 2
                          load(fullfile(cd, ['ChunkQS' num2str(pat(1,i)) '_win' num2str(win(1,j)) '_overlap' num2str(overlap(1,k)) '.mat']))  % load QS and AS; data and time vector of all AS/QS
                          
                           % creating time vector for QS
                           for l=1:length(RpeakQS)
                               if isempty(ChunkQS{1,l}) == 0 
                                  timingQS{1,l}=ChunkQS{1,l}(1,RpeakQS{1,l}(2:end)); % get the timestamp from AS time vector at the Rpeaks in s
                               end
                           end                          
                       end

                     

                        %lomb scargel
                        for l=1:length(RRdistanceAS)
                            if isempty(RRdistanceAS{1,l})==0
                                [pxxAS{1,l},fAS{1,l}] = plomb(RRdistanceAS{1,l},timingAS{1,l});
                                %convert to dB / frequency
                             %   pxxASdB=10*log10(pxxAS{1,l});
                            end
                        end
                        %normalizing
%                         maximum=max(pxx);
%                         pxx=pxx./maximum;
%                         pxx=pxx/length(pxx);
                        
                        
                        %plotting
                        if plotting  
                            for l=1:length(RRdistanceAS)
                                if isempty(RRdistanceAS{1,l})==0
                                    figure l
                                    plomb(RRdistanceAS{1,l},timingAS{1,l})                         
                                    title (['Lomb-Scargle periodogram AS patient: ' num2str(pat) ' win ' num2str(win) ' overlap ' num2str(overlap)])
        %                             xlabel('Frequency')
                                end
                            end
                        end
                   end
                                               
                   
                   
                   
%%%%%%%%%%%%%%%%QS
                    cd('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\windowed_ECG_chunks\Ralphs\')
                    if exist(fullfile(cd, ['RRdistanceQS' num2str(pat(1,i)) ' _ ' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,k)) '.mat']), 'file') == 2 % ==> 0 or 2
                        load(fullfile(cd, ['RRdistanceQS' num2str(pat(1,i)) ' _ ' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,k)) '.mat'])) % load  RR differnces     
                          
                                               
                        for l=1:length(RRdistanceQS)
                            RRdistanceQS{1,l}(any(isnan(RRdistanceQS{1,l})))=[]; %removing nans
                        end
                      
                        % removing outlier
                        for l=1:length(RRdistanceQS)
                            if isempty(RRdistanceQS{1,l})==0
                                all_idx = 1:length(RRdistanceQS{1,l});
                                outlier_idx = abs(RRdistanceQS{1,l} - median(RRdistanceQS{1,l})) > 3*std(RRdistanceQS{1,l}); % Find outlier idx
                                RRdistanceQS{1,l}(outlier_idx) = interp1(all_idx(~outlier_idx), RRdistanceQS{1,l}(~outlier_idx), all_idx(outlier_idx)); % Linearly interpolate over outlier idx for x
                            end
                        end
                        
                                                
                        %lomb scargel
                        for l=1:length(RRdistanceQS)
                            if isempty(RRdistanceQS{1,l})==0
                                [pxxQS{1,l},fQS{1,l}] = plomb(RRdistanceQS{1,l},timingQS{1,l});
                                %convert to dB / frequency
        %                         pxxQSdB=10*log10(pxxQS{1,l});
                            end
                        end
                        %normalizing
%                         maximum=max(pxx);
%                         pxx=pxx./maximum;
%                         pxx=pxx./length(pxx);
                        
                        %plotting
                        if plotting  
                            for l=1:length(RRdistanceQS)
                                if isempty(RRdistanceQS{1,l})==0
                                    figure l
                                    plomb(RRdistanceQS{1,l},timingQS{1,l})                         
                                    title (['Lomb-Scargle periodogram QS patient: ' num2str(pat) ' win ' num2str(win) ' overlap ' num2str(overlap)])
        %                             xlabel('Frequency')
                                end
                            end
                        end
                        
                                               
                    end
                    %saving
                    if saving
                        folder=('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\windowed_ECG_chunks\Ralphs\');
                        %saving power and frequency in mat file
                        save([folder 'PowerspectrumAS_' num2str(pat(i)) '_' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,k)) '.mat'],'pxxAS', 'fAS');
                        save([folder 'PowerspectrumQS_' num2str(pat(i)) '_' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,k)) '.mat'],'pxxQS', 'fQS');
                    end

                end
                    
            end
        end
    
   
    
end





