function pNNx(RR,Neonate,saving,savefolder,win,Session,S)
%Input
% RR: 5min RR distance data
% Neonate: Which patient
% saving: If saving is whished
% savefolder: Where to save
% win: Duration of the HRV window. Comon is 5min/300s

clearvars NN50 NN30 NN20 NN10 NN50 NN30 NN20 NN10 laenge


%%%%%%%%%%%%%%AS
pNN50(1:length(RR))=nan;pNN30=pNN50;pNN20=pNN50;pNN10=pNN50; %preallocation
for l=1:length(RR)
  laenge(l)=length(RR{l,1});
end
  laenge=sum(laenge);

  for i=1:length(RR)            

          pNN50(1,i)=(sum(abs(diff(RR{i,1}))>=50)/length(RR{i,1}))*100;
          if saving                     %saving R peaks positions in mat file                 
             Saving(pNN50,savefolder,Neonate,win,Session,S) 
          end% end if saving    

          pNN30(1,i)=(sum(abs(diff(RR{i,1}))>=30)/length(RR{i,1}))*100;
          if saving                     %saving R peaks positions in mat file                 
             Saving(pNN30,savefolder,Neonate,win,Session,S) 
          end% end if saving   

          pNN20(1,i)=(sum(abs(diff(RR{i,1}))>=20)/length(RR{i,1}))*100;
          if saving                     %saving R peaks positions in mat file                 
             Saving(pNN20,savefolder,Neonate,win,Session,S) 
          end% end if saving   

          pNN10(1,i)=(sum(abs(diff(RR{i,1}))>=10)/length(RR{i,1}))*100;
          if saving                     %saving R peaks positions in mat file                 
             Saving(pNN10,savefolder,Neonate,win,Session,S) 
          end% end if saving   
  end

      
    
        
  
%% %%%%%%%%%%replace 0 with 1337
%%%%%%%%%%% AS QS
%             if exist ('pNN50_AS','var')
%                  pNN50_AS(pNN50_AS==0)=1337; %all zeroes to 1337 to avoid confusion between AS and QS
%             end
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
 