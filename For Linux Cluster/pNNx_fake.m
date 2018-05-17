

function pNNx(pat,saving,win)
      for j=1:length(win)
         Neonate=pat;
         clearvars NN50_AS NN30_AS NN20_AS NN10_AS NN50_QS NN30_QS NN20_QS NN10_QS laenge

        %checking if file exist / loading
%         cd('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\ECG_30secepochs\Rpeaks\')
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
%%%%%%%%%%%%%%AS
        if exist('RRdistanceAS_win', 'var')==1
          for l=1:length(RRdistanceAS_win)
              laengeAS(l)=length(RRdistanceAS_win{1,l});
          end
              laengeAS=sum(laengeAS);

              for i=1:length(RRdistanceAS_win)            
                  if RRdistanceAS_win{1,i}
                      pNN50_AS(1,i)=sum(diff(RRdistanceAS_win{1,i})>=0.05)/laengeAS;
                      pNN30_AS(1,i)=sum(diff(RRdistanceAS_win{1,i})>=0.03 & diff(RRdistanceAS_win{1,i})<0.05)/laengeAS;
                      pNN20_AS(1,i)=sum(diff(RRdistanceAS_win{1,i})>=0.02 & diff(RRdistanceAS_win{1,i})<0.03)/laengeAS;
                      pNN10_AS(1,i)=sum(diff(RRdistanceAS_win{1,i})>=0.01 & diff(RRdistanceAS_win{1,i})<0.02)/laengeAS;
                  end
              end
        end
                          
%%%%%%%%%%%%%QS         
        if exist('RRdistanceQS_win','var')==1
            for l=1:length(RRdistanceQS_win)
              laengeQS(l)=length(RRdistanceQS_win{1,l});
          end
              laengeQS=sum(laengeQS);
              
           for i=1:length(RRdistanceQS_win)          
              if RRdistanceQS_win{1,i}
                  pNN50_QS(1,i)=(sum(diff(RRdistanceQS_win{1,i})>=0.05))/laengeQS;
                  pNN30_QS(1,i)=sum(diff(RRdistanceQS_win{1,i})>=0.03 & diff(RRdistanceQS_win{1,i})<0.05)/laengeQS;
                  pNN20_QS(1,i)=sum(diff(RRdistanceQS_win{1,i})>=0.02 & diff(RRdistanceQS_win{1,i})<0.03)/laengeQS;
                  pNN10_QS(1,i)=sum(diff(RRdistanceQS_win{1,i})>=0.01 & diff(RRdistanceQS_win{1,i})<0.02)/laengeQS;
              end
           end
        end
        
%%%%%%%%%%%%% Active Wake
        if exist('RRdistanceAalertness_win','var')==1
            for l=1:length(RRdistanceAalertness_win)
              laengeAalertness(l)=length(RRdistanceAalertness_win{1,l});
          end
              laengeAalertness=sum(laengeAalertness);
              
           for i=1:length(RRdistanceAalertness_win)          
              if RRdistanceAalertness_win{1,i}
                  pNN50_Aalertness(1,i)=(sum(diff(RRdistanceAalertness_win{1,i})>=0.05))/laengeAalertness;
                  pNN30_Aalertness(1,i)=sum(diff(RRdistanceAalertness_win{1,i})>=0.03 & diff(RRdistanceAalertness_win{1,i})<0.05)/laengeAalertness;
                  pNN20_Aalertness(1,i)=sum(diff(RRdistanceAalertness_win{1,i})>=0.02 & diff(RRdistanceAalertness_win{1,i})<0.03)/laengeAalertness;
                  pNN10_Aalertness(1,i)=sum(diff(RRdistanceAalertness_win{1,i})>=0.01 & diff(RRdistanceAalertness_win{1,i})<0.02)/laengeAalertness;
              end
           end
        end
%%%%%%%%%%%%% Quiet Wake
        if exist('RRdistanceQalertness_win','var')==1
            for l=1:length(RRdistanceQalertness_win)
              laengeQalertness(l)=length(RRdistanceQalertness_win{1,l});
          end
              laengeQalertness=sum(laengeQalertness);
              
           for i=1:length(RRdistanceQalertness_win)          
              if RRdistanceQalertness_win{1,i}
                  pNN50_Qalertness(1,i)=(sum(diff(RRdistanceQalertness_win{1,i})>=0.05))/laengeQalertness;
                  pNN30_Qalertness(1,i)=sum(diff(RRdistanceQalertness_win{1,i})>=0.03 & diff(RRdistanceQalertness_win{1,i})<0.05)/laengeQalertness;
                  pNN20_Qalertness(1,i)=sum(diff(RRdistanceQalertness_win{1,i})>=0.02 & diff(RRdistanceQalertness_win{1,i})<0.03)/laengeQalertness;
                  pNN10_Qalertness(1,i)=sum(diff(RRdistanceQalertness_win{1,i})>=0.01 & diff(RRdistanceQalertness_win{1,i})<0.02)/laengeQalertness;
              end
           end
        end   
        
%%%%%%%%%%%%% Transition
        if exist('RRdistancetransition_win','var')==1
            for l=1:length(RRdistancetransition_win)
              laengetransition(l)=length(RRdistancetransition_win{1,l});
          end
              laengetransition=sum(laengetransition);
              
           for i=1:length(RRdistancetransition_win)          
              if RRdistancetransition_win{1,i}
                  pNN50_transition(1,i)=(sum(diff(RRdistancetransition_win{1,i})>=0.05))/laengetransition;
                  pNN30_transition(1,i)=sum(diff(RRdistancetransition_win{1,i})>=0.03 & diff(RRdistancetransition_win{1,i})<0.05)/laengetransition;
                  pNN20_transition(1,i)=sum(diff(RRdistancetransition_win{1,i})>=0.02 & diff(RRdistancetransition_win{1,i})<0.03)/laengetransition;
                  pNN10_transition(1,i)=sum(diff(RRdistancetransition_win{1,i})>=0.01 & diff(RRdistancetransition_win{1,i})<0.02)/laengetransition;
              end
           end
        end  
        
%%%%%%%%%%%%% Position
        if exist('RRdistanceposition_win','var')==1
            for l=1:length(RRdistanceposition_win)
              laengeposition(l)=length(RRdistanceposition_win{1,l});
          end
              laengeposition=sum(laengeposition);
              
           for i=1:length(RRdistanceposition_win)          
              if RRdistanceposition_win{1,i}
                  pNN50_position(1,i)=(sum(diff(RRdistanceposition_win{1,i})>=0.05))/laengeposition;
                  pNN30_position(1,i)=sum(diff(RRdistanceposition_win{1,i})>=0.03 & diff(RRdistanceposition_win{1,i})<0.05)/laengeposition;
                  pNN20_position(1,i)=sum(diff(RRdistanceposition_win{1,i})>=0.02 & diff(RRdistanceposition_win{1,i})<0.03)/laengeposition;
                  pNN10_position(1,i)=sum(diff(RRdistanceposition_win{1,i})>=0.01 & diff(RRdistanceposition_win{1,i})<0.02)/laengeposition;
              end
           end
        end      
%%%%%%%%%%%%% Not reliable
        if exist('RRdistanceNot_reliable_win','var')==1
            for l=1:length(RRdistanceNot_reliable_win)
              laengeNot_reliable(l)=length(RRdistanceNot_reliable_win{1,l});
          end
              laengeNot_reliable=sum(laengeNot_reliable);
              
           for i=1:length(RRdistanceNot_reliable_win)          
              if RRdistanceNot_reliable_win{1,i}
                  pNN50_Not_reliable(1,i)=sum(diff(RRdistanceNot_reliable_win{1,i})>=0.05)/laengeNot_reliable;
                  pNN30_Not_reliable(1,i)=sum(diff(RRdistanceNot_reliable_win{1,i})>=0.03 & diff(RRdistanceNot_reliable_win{1,i})<0.05)/laengeNot_reliable;
                  pNN20_Not_reliable(1,i)=sum(diff(RRdistanceNot_reliable_win{1,i})>=0.02 & diff(RRdistanceNot_reliable_win{1,i})<0.03)/laengeNot_reliable;
                  pNN10_Not_reliable(1,i)=sum(diff(RRdistanceNot_reliable_win{1,i})>=0.01 & diff(RRdistanceNot_reliable_win{1,i})<0.02)/laengeNot_reliable;
              end
           end
        end   
%%%%%%%%%%%%replace 0 with 1337
%%%%%%%%%%% AS QS
            if exist ('pNN50_AS','var')
                 pNN50_AS(pNN50_AS==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end
            if exist ('pNN50_QS','var')
                 pNN50_QS(pNN50_QS==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end
            if exist ('pNN30_AS','var')
                 pNN30_AS(pNN30_AS==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end
            if exist ('pNN30_QS','var')
                 pNN30_QS(pNN30_QS==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end        
            if exist ('pNN20_AS','var')
                 pNN20_AS(pNN20_AS==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS 
            end
            if exist ('pNN20_QS','var')
                 pNN20_QS(pNN20_QS==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end    
            if exist ('pNN10_AS','var')
                 pNN10_AS(pNN10_AS==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end
            if exist ('pNN10_QS','var')
                 pNN10_QS(pNN10_QS==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end 
%%%%%%%%%%%% Wake
            if exist ('pNN50_Aalertness','var')
                 pNN50_Aalertness(pNN50_Aalertness==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end
            if exist ('pNN50_Qalertness','var')
                 pNN50_Qalertness(pNN50_Qalertness==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end
            if exist ('pNN30_Aalertness','var')
                 pNN30_Aalertness(pNN30_Aalertness==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end
            if exist ('pNN30_Qalertness','var')
                 pNN30_Qalertness(pNN30_Qalertness==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end        
            if exist ('pNN20_Aalertness','var')
                 pNN20_Aalertness(pNN20_Aalertness==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS 
            end
            if exist ('pNN20_Qalertness','var')
                 pNN20_Qalertness(pNN20_Qalertness==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end    
            if exist ('pNN10_Aalertness','var')
                 pNN10_Aalertness(pNN10_Aalertness==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end
            if exist ('pNN10_Qalertness','var')
                 pNN10_Qalertness(pNN10_Qalertness==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end 
%%%%%%%%% transition            
            if exist ('pNN50_transition','var')
                 pNN50_transition(pNN50_transition==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end

            if exist ('pNN30_transition','var')
                 pNN30_transition(pNN30_transition==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end

            if exist ('pNN20_transition','var')
                 pNN20_transition(pNN20_transition==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end    
            if exist ('pNN10_transition','var')
                 pNN10_transition(pNN10_transition==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end
%%%%%%%%% position           
            if exist ('pNN50_position','var')
                 pNN50_position(pNN50_position==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end

            if exist ('pNN30_position','var')
                 pNN30_position(pNN30_position==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end

            if exist ('pNN20_position','var')
                 pNN20_position(pNN20_position==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end    
            if exist ('pNN10_position','var')
                 pNN10_position(pNN10_position==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end          

%%%%%%%%% Not reliable           
            if exist ('pNN50_Not_reliable','var')
                 pNN50_Not_reliable(pNN50_Not_reliable==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end

            if exist ('pNN30_Not_reliable','var')
                 pNN30_Not_reliable(pNN30_Not_reliable==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end

            if exist ('pNN20_Not_reliable','var')
                 pNN20_Not_reliable(pNN20_Not_reliable==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end    
            if exist ('pNN10_Not_reliable','var')
                 pNN10_Not_reliable(pNN10_Not_reliable==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
            end                      
%%%%%%%%%%%Saving

        if saving
            folder='C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\HRV features\';
%%%%%%%%% AS QS            
            if exist('pNN50_AS','var') & exist ('pNN50_QS','var') & exist('pNN30_AS','var') & exist ('pNN30_QS','var') & exist('pNN20_AS','var') & exist ('pNN20_QS','var') & exist('pNN10_AS','var') & exist ('pNN10_QS','var')
             save([folder 'pNNx_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'pNN50_AS','pNN30_AS','pNN20_AS','pNN10_AS','pNN50_QS', 'pNN30_QS','pNN20_QS', 'pNN10_QS');
            elseif exist('pNN50_AS')  & exist('pNN30_AS') &  exist('pNN20_AS') &  exist('pNN10_AS') 
             save([folder 'pNNx_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'pNN50_AS','pNN30_AS','pNN20_AS','pNN10_AS');
            elseif exist ('pNN50_QS') &  exist ('pNN30_QS') &  exist ('pNN20_QS') & exist ('pNN10_QS')
             save([folder 'pNNx_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'pNN50_QS', 'pNN30_QS','pNN20_QS', 'pNN10_QS');
            end
%%%%%%%% Wake            
            if exist('pNN50_Aalertness') & exist ('pNN50_Qalertness') & exist('pNN30_Aalertness') & exist ('pNN30_Qalertness') & exist('pNN20_Aalertness') & exist ('pNN20_Qalertness') & exist('pNN10_Aalertness') & exist ('pNN10_Qalertness')
             save([folder 'pNNx_Wake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'pNN50_Aalertness','pNN30_Aalertness','pNN20_Aalertness','pNN10_Aalertness','pNN50_Qalertness', 'pNN30_Qalertness','pNN20_Qalertness', 'pNN10_Qalertness');
            elseif exist('pNN50_Aalertness')  & exist('pNN30_Aalertness') &  exist('pNN20_Aalertness') &  exist('pNN10_Aalertness') 
             save([folder 'pNNx_Wake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'pNN50_Aalertness','pNN30_Aalertness','pNN20_Aalertness','pNN10_Aalertness');
            elseif exist ('pNN50_Qalertness') &  exist ('pNN30_Qalertness') &  exist ('pNN20_Qalertness') & exist ('pNN10_Qalertness')
             save([folder 'pNNx_Wake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'pNN50_Qalertness', 'pNN30_Qalertness','pNN20_Qalertness', 'pNN10_Qalertness');
            end                            
%%%%%%%% transition            
            if exist('pNN50_transition','var') & exist('pNN30_transition','var') & exist('pNN20_transition','var') & exist('pNN10_transition','var')
             save([folder 'pNNx_transition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'pNN50_transition','pNN30_transition','pNN20_transition','pNN10_transition','pNN50_transition', 'pNN30_transition','pNN20_transition', 'pNN10_transition');         
            end   
%%%%%%%% position            
            if exist('pNN50_position','var') & exist('pNN30_position','var') & exist('pNN20_position','var') & exist('pNN10_position','var')
             save([folder 'pNNx_position_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'pNN50_position','pNN30_position','pNN20_position','pNN10_position','pNN50_position', 'pNN30_position','pNN20_position', 'pNN10_position');         
            end  
%%%%%%%% Not reliable            
            if exist('pNN50_Not_reliable','var') & exist('pNN30_Not_reliable','var') & exist('pNN20_Not_reliable','var') & exist('pNN10_Not_reliable','var')
             save([folder 'pNNx_Not_reliable_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'pNN50_Not_reliable','pNN30_Not_reliable','pNN20_Not_reliable','pNN10_Not_reliable','pNN50_Not_reliable', 'pNN30_Not_reliable','pNN20_Not_reliable', 'pNN10_Not_reliable');         
            end   
        end% saving
        
      end %length
    end
