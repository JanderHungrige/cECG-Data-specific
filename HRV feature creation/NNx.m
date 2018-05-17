function NNx(RR,Neonate,saving,savefolder,win,Session,S)
%Input
% RR: 5min RR distance data
% Neonate: Which patient
% saving: If saving is whished
% savefolder: Where to save
% win: Duration of the HRV window. Comon is 5min/300s

clearvars NN50 NN30 NN20 NN10 NN50 NN30 NN20 NN10 
            
%%%%%%%%%%%%AS                        

NN50(1:length(RR))=nan;NN30=NN50;NN20=NN50;NN10=NN50; %preallocation
    for i=1:length(RR)
           
          NN50(1,i)=sum(abs(diff(RR{i,1}))>=50);
          if saving                     %saving R peaks positions in mat file                 
             Saving(NN50,savefolder,Neonate,win,Session,S) 
          end% end if saving  

          NN30(1,i)=sum(abs(diff(RR{i,1}))>=30);
          if saving                     %saving R peaks positions in mat file                 
             Saving(NN30,savefolder,Neonate,win,Session,S) 
          end% end if saving

          NN20(1,i)=sum(abs(diff(RR{i,1}))>=20);
          if saving                     %saving R peaks positions in mat file                 
             Saving(NN20,savefolder,Neonate,win,Session,S) 
          end% end if saving   

          NN10(1,i)=sum(abs(diff(RR{i,1}))>=10);
          if saving                     %saving R peaks positions in mat file                 
             Saving(NN10,savefolder,Neonate,win,Session,S) 
          end% end if saving 
          
    end

      
    
                   
                    
%% %%%%%%%%%%%%replace 0 with 1337
% %%%%%%%%%% AS QS
%             if exist ('NN50_AS','var')
%                  NN50_AS(NN50_AS==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
%             end
%    
           

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

