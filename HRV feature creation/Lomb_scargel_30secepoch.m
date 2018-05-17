%Frequency analysis

%[pxx,f] = plomb(x,t) returns the Lomb-Scargle power spectral density (PSD) estimate, pxx, of a signal, x, that is sampled at the instants specified in t. 
% t must increase monotonically but need not be uniformly spaced. All elements of t must be nonnegative. pxx is evaluated at the frequencies returned in f.


function Lomb_scargel_30secepoch(Neonate,plotting,win,saving)
%     for i=1:length(Neonate)
        for j=1:length(win)
            
%%%%%%%%Loading files
%%%%%% AS QS, Wake, Transition, Position, Not Reliable
        cd('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\R peak detection and HRV')
            k=2; %loading windowed ECG/RR data
            loading_RR_data(k,Neonate,win(1,j),win) %this function loads all the data
        folder=('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\R peak detection and HRV\');
            load(fullfile(folder,['loading_RR_data workspace']));   % The workspace has to be loaded from "loading_RR_data" otherwise the files do not appear in this workspace
            delete(fullfile(folder,['loading_RR_data workspace'])); % delete it after use
            
        
                       
%%%%%%%%Data preprocessing
                     if exist('RRdistanceAS_win','var') == 1
                %removing outlier
                        for l=1:length(RRdistanceAS_win)
                            if isempty(RRdistanceAS_win{1,l})==0
                                all_idx = 1:length(RRdistanceAS_win{1,l});
                                outlier_idx = abs(RRdistanceAS_win{1,l} - median(RRdistanceAS_win{1,l})) > 3*std(RRdistanceAS_win{1,l}); % Find outlier idx
                                RRdistanceAS_win{1,l}(outlier_idx) = interp1(all_idx(~outlier_idx), RRdistanceAS_win{1,l}(~outlier_idx), all_idx(outlier_idx)); % Linearly interpolate over outlier idx for x
                            end
                        end
                     end
                        
                       if exist('RRdistanceQS_win','var') == 1
                        % removing outlier
                        for l=1:length(RRdistanceQS_win)
                            if isempty(RRdistanceQS_win{1,l})==0
                                all_idx = 1:length(RRdistanceQS_win{1,l});
                                outlier_idx = abs(RRdistanceQS_win{1,l} - median(RRdistanceQS_win{1,l})) > 3*std(RRdistanceQS_win{1,l}); % Find outlier idx
                                RRdistanceQS_win{1,l}(outlier_idx) = interp1(all_idx(~outlier_idx), RRdistanceQS_win{1,l}(~outlier_idx), all_idx(outlier_idx)); % Linearly interpolate over outlier idx for x
                            end
                         end   
                       end
                           
                        
%%%%%%%%Creating vectors

%                 %creating time vector for AS
%                            for l=1:length(RpeakAS_win)
%                                if isempty(ChunkAS_win{1,l}) == 0 
%                                   timingAS{1,l}=ChunkAS_win{1,l}(1,RpeakAS_win{1,l}(2:end)); % get the timestamp from AS time vector at the Rpeaks in s
%                                end
%                            end
%                        
%                           
%                 %creating time vector for QS
%                            for l=1:length(RpeakQS_win)
%                                if isempty(ChunkQS_win{1,l}) == 0 
%                                   timingQS{1,l}=ChunkQS_win{1,l}(1,RpeakQS_win{1,l}(2:end)); % get the timestamp from AS time vector at the Rpeaks in s
%                                end
%                            end                          
%                            
%%%%%%%%Create spectrum

                %lomb scargel
                %%%%%%%%%%%%% AS QS 
                    if exist('RRdistanceAS_win','var') == 1
                        for l=1:length(RRdistanceAS_win)
                            if isempty(RRdistanceAS_win{1,l})==0
                                [pxxAS{1,l},fAS{1,l}] = plomb(RRdistanceAS_win{1,l},RpeakAStiming_win{1,l});
                                %convert to dB / frequency
                                %pxxASdB=10*log10(pxxAS{1,l});
                            end
                        end
                    end
                    
                     if exist('RRdistanceQS_win','var') == 1                        
                        for l=1:length(RRdistanceQS_win)
                            if isempty(RRdistanceQS_win{1,l})==0
                                [pxxQS{1,l},fQS{1,l}] = plomb(RRdistanceQS_win{1,l},RpeakQStiming_win{1,l});
                            end
                        end  
                     end
                %%%%%%%%%%%%% WAKE
                     if exist('RRdistanceAalertness_win','var') == 1                                        
                        for l=1:length(RRdistanceAalertness_win)
                            if isempty(RRdistanceAalertness_win{1,l})==0
                                [pxxAalertness{1,l},fAalertness{1,l}] = plomb(RRdistanceAalertness_win{1,l},RpeakAalertnesstiming_win{1,l});
                            end
                        end
                     end
                     
                     if exist('RRdistanceQalertness_win','var') == 1  
                        for l=1:length(RRdistanceQalertness_win)
                            if isempty(RRdistanceQalertness_win{1,l})==0
                                [pxxQalertness{1,l},fQalertness{1,l}] = plomb(RRdistanceQalertness_win{1,l},RpeakQalertnesstiming_win{1,l});
                            end
                        end         
                     end
                %%%%%%%%%%%%% TRANSITION
                     if exist('RRdistancetransition_win','var') == 1                
                        for l=1:length(RRdistancetransition_win)
                            if isempty(RRdistancetransition_win{1,l})==0
                                [pxxtransition{1,l},ftransition{1,l}] = plomb(RRdistancetransition_win{1,l},Rpeaktransitiontiming_win{1,l});
                            end
                        end
                     end
                %%%%%%%%%%%%% POSITION 
                     if exist('RRdistanceposition_win','var') == 1                
                        for l=1:length(RRdistanceposition_win)
                            if isempty(RRdistanceposition_win{1,l})==0
                                [pxxposition{1,l},fposition{1,l}] = plomb(RRdistanceposition_win{1,l},Rpeakpositiontiming_win{1,l});
                            end
                        end   
                     end
                %%%%%%%%%%%%% NOT RELIABLE 
                     if exist('RRdistanceNot_reliable_win','var') == 1                
                        for l=1:length(RRdistanceNot_reliable_win)
                            if isempty(RRdistanceNot_reliable_win{1,l})==0
                                [pxxNot_reliable{1,l},fNot_reliable{1,l}] = plomb(RRdistanceNot_reliable_win{1,l},RpeakNot_reliabletiming_win{1,l});
                            end
                        end 
                     end

                                            
%%%%%%%%% saving
                    if saving
                        folder=('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\ECG_30secepochs\spectrum\');
 %%%%%%%%% AS QS                      
                       if exist('pxxAS','var') ==1 & exist('pxxQS','var') ==1 & exist('fAS','var') ==1 & exist('fQS','var') ==1
                           save([folder 'Powerspectrum30_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'pxxAS','pxxQS', 'fAS', 'fQS')  ; 
                       elseif exist('pxxAS','var') ==1 & exist('fAS','var') ==1
                           save([folder 'Powerspectrum30_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'pxxAS', 'fAS')  ; 
                       elseif exist('pxxQS','var') ==1 & exist('fQS','var') ==1
                           save([folder 'Powerspectrum30_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'pxxQS', 'fQS')  ; 
                       end  
 %%%%%%%%% WAKE                      
                       if exist('pxxAalertness','var') ==1 & exist('pxxQalertness','var') ==1 & exist('fAalertness','var') ==1 & exist('fQalertness','var') ==1
                           save([folder 'PowerspectrumAalertness_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'pxxAalertness','pxxQalertness', 'fAalertness', 'fQalertness')  ; 
                       elseif exist('pxxAalertness','var') ==1 & exist('fAalertness','var') ==1
                           save([folder 'PowerspectrumAalertness_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'pxxAalertness', 'fAalertness')  ; 
                       elseif exist('pxxQalertness','var') ==1 & exist('fQalertness','var') ==1
                           save([folder 'PowerspectrumQalertness_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'pxxQalertness', 'fQalertness')  ; 
                       end                         
%%%%%%%%% TRANSITION
                       if exist('pxxtransition','var') ==1 & exist('ftransition','var') ==1
                           save([folder 'PowerspectrumTransition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'pxxtransition', 'ftransition')  ; 
                       end 
%%%%%%%%% POSITION
                       if exist('pxxposition','var') ==1 & exist('fposition','var') ==1
                           save([folder 'PowerspectrumPosition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'pxxposition', 'fposition')  ; 
                       end 
%%%%%%%%% NOT RELIABLE
                       if exist('pxxNot_reliable','var') ==1 & exist('fNot_reliable','var') ==1
                           save([folder 'PowerspectrumNot_reliable_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'pxxNot_reliable', 'fNot_reliable')  ; 
                       end 
                       
                       
                    end  %if saving
                    
          
        end%win
%    end% Neonate
    
%%%%%%%% PLOTING
                                                
                        if plotting  
                            for l=1:length(RRdistanceAS_win)
                                if isempty(RRdistanceAS_win{1,l})==0
                                    figure l
                                    plomb(RRdistanceAS_win{1,l},RpeakAStiming_win{1,l})                         
                                    title (['Lomb-Scargle periodogram AS patient: ' num2str(Neonate) ' win ' num2str(win) ])
        %                             xlabel('Frequency')
                                end
                            end
                        end
                     
                        
                        if plotting  
                            for l=1:length(RRdistanceQS_win)
                                if isempty(RRdistanceQS_win{1,l})==0
                                    figure l
                                    plomb(RRdistanceQS_win{1,l},RpeakQStiming_win{1,l})                         
                                    title (['Lomb-Scargle periodogram QS patient: ' num2str(Neonate) ' win ' num2str(win(1,j))])
        %                             xlabel('Frequency')
                                end
                            end
                        end                                             
    
end





