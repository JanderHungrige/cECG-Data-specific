

function NNx(pat,saving,win)
 Neonate=pat;

       for j=1:length(win)
              clearvars NN50_AS NN30_AS NN20_AS NN10_AS NN50_QS NN30_QS NN20_QS NN10_QS

          %checking if file exist / loading
%          cd('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\ECG_30secepochs\Rpeaks\')
         cd('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\ECG_30secepochs\ECG clustered by annotations\unsorted')        
         if exist(fullfile(cd, ['RRdistance_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']), 'file') == 2 % ==> 0 or 2
            load(fullfile(cd, ['RRdistance_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']))
         end
        cd('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\ECG_30secepochs\ECG clustered by annotations\unsorted')        
         if exist(fullfile(cd, ['RRdistanceWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']), 'file') == 2 % ==> 0 or 2
            load(fullfile(cd, ['RRdistanceWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']))
         end       

        cd('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\ECG_30secepochs\ECG clustered by annotations\unsorted')        
         if exist(fullfile(cd, ['RRdistancetransition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']), 'file') == 2 % ==> 0 or 2
            load(fullfile(cd, ['RRdistancetransition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']))
         end 

        cd('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\ECG_30secepochs\ECG clustered by annotations\unsorted')        
         if exist(fullfile(cd, ['RRdistanceposition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']), 'file') == 2 % ==> 0 or 2
            load(fullfile(cd, ['RRdistanceposition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']))
         end   

        cd('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\ECG_30secepochs\ECG clustered by annotations\unsorted')        
         if exist(fullfile(cd, ['RRdistanceNot_reliable_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']), 'file') == 2 % ==> 0 or 2
            load(fullfile(cd, ['RRdistanceNot_reliable_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']))
        end       
            
%%%%%%%%%%%%AS                        
%                   for i=1:length(RRdistanceAS)
%                          RRdistanceAS{1,i}(any(isnan( RRdistanceAS{1,i})))=[]; %removing nans
%                       if RRdistanceAS{1,i}
%                           
%                           NN50_AS(1,i)=sum(diff(RRdistanceAS{1,i})>=0.05);
%                           NN30_AS(1,i)=sum(diff(RRdistanceAS{1,i})>=0.03);
%                           NN20_AS(1,i)=sum(diff(RRdistanceAS{1,i})>=0.02));
%                           NN10_AS(1,i)=sum(diff(RRdistanceAS{1,i})>=0.01);
%                       end
%                    end
            if exist('RRdistanceAS_win','var')==1
                    for i=1:length(RRdistanceAS_win)
                       if RRdistanceAS_win{1,i}
                          NN50_AS(1,i)=sum(diff(RRdistanceAS_win{1,i})>=0.05);
                          NN30_AS(1,i)=sum(diff(RRdistanceAS_win{1,i})>=0.03 & diff(RRdistanceAS_win{1,i}) <0.05);
                          NN20_AS(1,i)=sum(diff(RRdistanceAS_win{1,i})>=0.02 & diff(RRdistanceAS_win{1,i}) <0.03);
                          NN10_AS(1,i)=sum(diff(RRdistanceAS_win{1,i})>=0.01 & diff(RRdistanceAS_win{1,i}) <0.02);
                      end
                    end
            end
%%%%%%%%%%%%%QS
            if exist('RRdistanceQS_win','var')==1
                   for i=1:length(RRdistanceQS_win)
                      if RRdistanceQS_win{1,i}
                          NN50_QS(1,i)=sum(diff(RRdistanceQS_win{1,i})>=0.05);
                          NN30_QS(1,i)=sum(diff(RRdistanceQS_win{1,i})>=0.03 & diff(RRdistanceQS_win{1,i})<0.05);
                          NN20_QS(1,i)=sum(diff(RRdistanceQS_win{1,i})>=0.02 & diff(RRdistanceQS_win{1,i})<0.03);
                          NN10_QS(1,i)=sum(diff(RRdistanceQS_win{1,i})>=0.01 & diff(RRdistanceQS_win{1,i})<0.02);
                      end
                   end
            end
%%%%%%%%%%%%% WAKE
            if exist('RRdistanceAalertness_win','var')==1
                   for i=1:length(RRdistanceAalertness_win)
                      if RRdistanceAalertness_win{1,i}
                          NN50_Aalertness(1,i)=sum(diff(RRdistanceAalertness_win{1,i})>=0.05);
                          NN30_Aalertness(1,i)=sum(diff(RRdistanceAalertness_win{1,i})>=0.03 & diff(RRdistanceAalertness_win{1,i})<0.05);
                          NN20_Aalertness(1,i)=sum(diff(RRdistanceAalertness_win{1,i})>=0.02 & diff(RRdistanceAalertness_win{1,i})<0.03);
                          NN10_Aalertness(1,i)=sum(diff(RRdistanceAalertness_win{1,i})>=0.01 & diff(RRdistanceAalertness_win{1,i})<0.02);
                      end
                   end
            end            
            if exist('RRdistanceQalertness_win','var')==1
                   for i=1:length(RRdistanceQalertness_win)
                      if RRdistanceQalertness_win{1,i}
                          NN50_Qalertness(1,i)=sum(diff(RRdistanceQalertness_win{1,i})>=0.05);
                          NN30_Qalertness(1,i)=sum(diff(RRdistanceQalertness_win{1,i})>=0.03 & diff(RRdistanceQalertness_win{1,i})<0.05);
                          NN20_Qalertness(1,i)=sum(diff(RRdistanceQalertness_win{1,i})>=0.02 & diff(RRdistanceQalertness_win{1,i})<0.03);
                          NN10_Qalertness(1,i)=sum(diff(RRdistanceQalertness_win{1,i})>=0.01 & diff(RRdistanceQalertness_win{1,i})<0.02);
                      end
                   end
            end      
%%%%%%%%%%%%% TRANSITION
            if exist('RRdistancetransition_win','var')==1
                   for i=1:length(RRdistancetransition_win)
                      if RRdistancetransition_win{1,i}
                          NN50_transition(1,i)=sum(diff(RRdistancetransition_win{1,i})>=0.05);
                          NN30_transition(1,i)=sum(diff(RRdistancetransition_win{1,i})>=0.03 & diff(RRdistancetransition_win{1,i})<0.05);
                          NN20_transition(1,i)=sum(diff(RRdistancetransition_win{1,i})>=0.02 & diff(RRdistancetransition_win{1,i})<0.03);
                          NN10_transition(1,i)=sum(diff(RRdistancetransition_win{1,i})>=0.01 & diff(RRdistancetransition_win{1,i})<0.02);
                      end
                   end
            end        
%%%%%%%%%%%%% POSITION
            if exist('RRdistanceposition_win','var')==1
                   for i=1:length(RRdistanceposition_win)
                      if RRdistanceposition_win{1,i}
                          NN50_position(1,i)=sum(diff(RRdistanceposition_win{1,i})>=0.05);
                          NN30_position(1,i)=sum(diff(RRdistanceposition_win{1,i})>=0.03 & diff(RRdistanceposition_win{1,i})<0.05);
                          NN20_position(1,i)=sum(diff(RRdistanceposition_win{1,i})>=0.02 & diff(RRdistanceposition_win{1,i})<0.03);
                          NN10_position(1,i)=sum(diff(RRdistanceposition_win{1,i})>=0.01 & diff(RRdistanceposition_win{1,i})<0.02);
                      end
                   end
%%%%%%%%%%%%% NOT RELIABLE
            if exist('RRdistanceNot_reliable_win','var')==1
                   for i=1:length(RRdistanceNot_reliable_win)
                      if RRdistanceNot_reliable_win{1,i}
                          NN50_Not_reliable(1,i)=sum(diff(RRdistanceNot_reliable_win{1,i})>=0.05);
                          NN30_Not_reliable(1,i)=sum(diff(RRdistanceNot_reliable_win{1,i})>=0.03 & diff(RRdistanceNot_reliable_win{1,i})<0.05);
                          NN20_Not_reliable(1,i)=sum(diff(RRdistanceNot_reliable_win{1,i})>=0.02 & diff(RRdistanceNot_reliable_win{1,i})<0.03);
                          NN10_Not_reliable(1,i)=sum(diff(RRdistanceNot_reliable_win{1,i})>=0.01 & diff(RRdistanceNot_reliable_win{1,i})<0.02);
                      end
                   end                   
            end             
%%%%%%%%%%%%replace 0 with 1337
%%%%%%%%%% AS QS
            if exist ('NN50_AS','var')
                 NN50_AS(NN50_AS==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end
            if exist ('NN50_QS','var')
                 NN50_QS(NN50_QS==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end
            if exist ('NN30_AS','var')
                 NN30_AS(NN30_AS==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end
            if exist ('NN30_QS','var')
                 NN30_QS(NN30_QS==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end        
            if exist ('NN20_AS','var')
                 NN20_AS(NN20_AS==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS 
            end
            if exist ('NN20_QS','var')
                 NN20_QS(NN20_QS==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end    
            if exist ('NN10_AS','var')
                 NN10_AS(NN10_AS==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end
            if exist ('NN10_QS','var')
                 NN10_QS(NN10_QS==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end     
%%%%%%%%%%%% Wake
            if exist ('NN50_Aalertness','var')
                 NN50_Aalertness(NN50_Aalertness==0)=1337; %all zeroes to nan to avoid confusion between AS and QS
            end
            if exist ('NN50_Qalertness','var')
                 NN50_Qalertness(NN50_Qalertness==0)=1337; %all zeroes to nan to avoid confusion between AS and QS
            end
            if exist ('NN30_Aalertness','var')
                 NN30_Aalertness(NN30_Aalertness==0)=1337; %all zeroes to nan to avoid confusion between AS and QS
            end
            if exist ('pNN30_Qalertness','var')
                 NN30_Qalertness(NN30_Qalertness==0)=1337; %all zeroes to nan to avoid confusion between AS and QS
            end        
            if exist ('pNN20_Aalertness','var')
                 NN20_Aalertness(NN20_Aalertness==0)=1337; %all zeroes to nan to avoid confusion between AS and QS 
            end
            if exist ('pNN20_Qalertness','var')
                 NN20_Qalertness(NN20_Qalertness==0)=1337; %all zeroes to nan to avoid confusion between AS and QS
            end    
            if exist ('pNN10_Aalertness','var')
                 NN10_Aalertness(NN10_Aalertness==0)=1337; %all zeroes to nan to avoid confusion between AS and QS
            end
            if exist ('pNN10_Qalertness','var')
                 NN10_Qalertness(NN10_Qalertness==0)=1337; %all zeroes to nan to avoid confusion between AS and QS
            end 
%%%%%%%%% transition            
            if exist ('pNN50_transition','var')
                 NN50_transition(NN50_transition==0)=1337; %all zeroes to nan to avoid confusion between AS and QS
            end

            if exist ('NN30_transition','var')
                 NN30_transition(NN30_transition==0)=1337; %all zeroes to nan to avoid confusion between AS and QS
            end

            if exist ('NN20_transition','var')
                 NN20_transition(NN20_transition==0)=1337; %all zeroes to nan to avoid confusion between AS and QS
            end    
            if exist ('NN10_transition','var')
                 NN10_transition(NN10_transition==0)=1337; %all zeroes to nan to avoid confusion between AS and QS
            end
%%%%%%%%% position           
            if exist ('NN50_position','var')
                 NN50_position(NN50_position==0)=1337; %all zeroes to nan to avoid confusion between AS and QS
            end

            if exist ('NN30_position','var')
                 NN30_position(NN30_position==0)=1337; %all zeroes to nan to avoid confusion between AS and QS
            end

            if exist ('NN20_position','var')
                 NN20_position(NN20_position==0)=1337; %all zeroes to nan to avoid confusion between AS and QS
            end    
            if exist ('NN10_position','var')
                 NN10_position(NN10_position==0)=1337; %all zeroes to nan to avoid confusion between AS and QS
            end          

%%%%%%%%% Not reliable           
            if exist ('NN50_Not_reliable','var')
                 NN50_Not_reliable(NN50_Not_reliable==0)=1337; %all zeroes to nan to avoid confusion between AS and QS
            end

            if exist ('NN30_Not_reliable','var')
                 NN30_Not_reliable(NN30_Not_reliable==0)=1337; %all zeroes to nan to avoid confusion between AS and QS
            end

            if exist ('NN20_Not_reliable','var')
                 NN20_Not_reliable(NN20_Not_reliable==0)=1337; %all zeroes to nan to avoid confusion between AS and QS
            end    
            if exist ('NN10_Not_reliable','var')
                 NN10_Not_reliable(NN10_Not_reliable==0)=1337; %all zeroes to nan to avoid confusion between AS and QS
            end                            
%%%%%%%%%%%%Saving
            if saving
                folder='C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\HRV features\';
%%%%%%%% AS QS                
                if exist('NN50_AS') & exist ('NN50_QS')& exist('NN30_AS') & exist ('NN30_QS') & exist('NN20_AS') & exist ('NN20_QS') & exist('NN10_AS') & exist ('NN10_QS')
                 save([folder 'NNx_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'NN50_AS','NN30_AS','NN20_AS','NN10_AS','NN50_QS', 'NN30_QS','NN20_QS', 'NN10_QS');
                elseif  exist('NN50_AS')  & exist('NN30_AS')  & exist('NN20_AS') & exist('NN10_AS')
                 save([folder 'NNx_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'NN50_AS','NN30_AS','NN20_AS','NN10_AS');
                elseif  exist ('NN50_QS') &  exist ('NN30_QS') &  exist ('NN20_QS') &  exist ('NN10_QS')
                 save([folder 'NNx_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'NN50_QS', 'NN30_QS','NN20_QS', 'NN10_QS');
                end
%%%%%%%% Wake            
            if exist('NN50_Aalertness') & exist ('NN50_Qalertness') & exist('NN30_Aalertness') & exist ('NN30_Qalertness') & exist('NN20_Aalertness') & exist ('NN20_Qalertness') & exist('NN10_Aalertness') & exist ('NN10_Qalertness')
             save([folder 'NNx_Wake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'NN50_Aalertness','NN30_Aalertness','NN20_Aalertness','NN10_Aalertness','NN50_Qalertness', 'NN30_Qalertness','NN20_Qalertness', 'NN10_Qalertness');
            elseif exist('NN50_Aalertness')  & exist('NN30_Aalertness') &  exist('NN20_Aalertness') &  exist('NN10_Aalertness') 
             save([folder 'NNx_Wake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'NN50_Aalertness','NN30_Aalertness','NN20_Aalertness','NN10_Aalertness');
            elseif exist ('NN50_Qalertness') &  exist ('NN30_Qalertness') &  exist ('NN20_Qalertness') & exist ('NN10_Qalertness')
             save([folder 'NNx_Wake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'NN50_Qalertness', 'NN30_Qalertness','NN20_Qalertness', 'NN10_Qalertness');
            end                            
%%%%%%%% transition            
            if exist('NN50_transition','var') & exist('NN30_transition','var') & exist('NN20_transition','var') & exist('NN10_transition','var')
             save([folder 'NNx_transition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'NN50_transition','NN30_transition','NN20_transition','NN10_transition','NN50_transition', 'NN30_transition','NN20_transition', 'NN10_transition');         
            end   
%%%%%%%% position            
            if exist('NN50_position','var') & exist('NN30_position','var') & exist('NN20_position','var') & exist('NN10_position','var')
             save([folder 'NNx_position_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'NN50_position','NN30_position','NN20_position','NN10_position','NN50_position', 'NN30_position','NN20_position', 'NN10_position');         
            end  
%%%%%%%% Not reliable            
            if exist('NN50_Not_reliable','var') & exist('NN30_Not_reliable','var') & exist('NN20_Not_reliable','var') & exist('NN10_Not_reliable','var')
             save([folder 'NNx_Not_reliable_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'NN50_Not_reliable','NN30_Not_reliable','NN20_Not_reliable','NN10_Not_reliable','NN50_Not_reliable', 'NN30_Not_reliable','NN20_Not_reliable', 'NN10_Not_reliable');         
            end   
                
       
            end   % saving               
      
      end%win
end