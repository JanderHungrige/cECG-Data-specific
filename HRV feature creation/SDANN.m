function SDANN(RR,Neonate,saving,savefolder,win,faktor,S) 
%Input
% RR: 5min RR distance data
% Neonate: Which patient
% saving: If saving is whished
% savefolder: Where to save
% win: Duration of the HRV window. Comon is 5min/300s
% factor: How much the SIgnal is shifted each itteration. 30s is common. 

% #1 first calculate STD for each 30s Epoch.
% #2 merge those 30s Epoch to 5min
% #3 mean them 
 
% #1
for i=1:length(RR)
  if RR{i,1}~=0
  SDNN_30{1,i}=nanstd(RR{i,1});
  end
end     
 
% #2
win_jumps=faktor/faktor; % as it is already on 30s epoch we need a shift by 1
Fenster=win/30; % 300/30= 10 parts a 30s => 5min 

m=1;
uebrig=length(SDNN_30);  % how many minutes are left           
for k=1:win_jumps:length(SDNN_30)
   if k+Fenster<length(SDNN_30) 
       SDANN{1,m}=SDNN_30(1,k:k+Fenster-1); 
   elseif k+Fenster>=length(SDNN_30) && win_jumps<=uebrig && k>win_jumps*(Fenster-(uebrig/win_jumps)*win_jumps)/win_jumps
       rechts=uebrig/win_jumps;% How many epochs are still left 
       links=(Fenster-rechts*win_jumps)/win_jumps; % how many epochs do we have to atahe from the left to get a full 300s window
       SDANN{1,m}=SDNN_30(1,k-win_jumps*links:k+win_jumps*rechts-1);
   elseif k+Fenster>=length(SDNN_30) &&  win_jumps>uebrig && k>win_jumps*(Fenster-(uebrig/win_jumps)*win_jumps)/win_jumps     
       rechts=uebrig/win_jumps;% How many epochs are still left 
       links=(Fenster-rechts*win_jumps)/win_jumps; % how many epochs do we have to atahe from the left to get a full 300s window
       SDANN{1,m}=SDNN_30(1,k-win_jumps*links:end);       
   else
       SDANN{1,m}=SDNN_30(1,k:end);
%            break       % if you want to end with the same length for the clast cell elementas the others use break. But than the ECG_win_300 is one element shorter thatn ECG_win_30    
   end
   uebrig=length(SDNN_30)-(k+win_jumps);  % how many minutes are left           
   m=m+1;
end

% Original &&&&&&&&&&&&&&&&&&&&&&&&& Working
% for k=1:win_jumps:length(SDNN_30)
%    if k+Fenster<length(SDNN_30) 
%      SDANN{1,M}=SDNN_30(1,k:k+Fenster-1);      
%    elseif k+Fenster>=length(SDNN_30) 
%        SDANN{1,M}=SDNN_30(1,k:end);
%        break
%    end
%    M=M+1;
% end
% Original &&&&&&&&&&&&&&&&&&&&&&&&& Working

% #3   
for N=1:length(SDANN)
   SDANN{1,N}=nanmean(cell2mat(SDANN{1,N}));
end

            
%%%%%%%%%%%% SAVING            
    if saving                     %saving R peaks positions in mat file                 
       Saving(SDANN,savefolder,Neonate,win,S) 
    end% end if saving 
    
    
  
end

%% Nested saving
    function Saving(Feature,savefolder, Neonate,win,S)
        if exist('Feature','var')==1
            name=inputname(1); % variable name of function input
            save([savefolder name '_Session_' num2str(S) 'pat_' num2str(Neonate)],'Feature')
        else
            disp(['saving of ' name ' not possible'])
        end       
    end
 
 
 
