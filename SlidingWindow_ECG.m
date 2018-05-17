function [ECG_win_300,ECG_win_30,t_300,t_30]=SlidingWindow_ECG(ECG,t_ECG,Neonate,saving,folder,factor,win,S)  
% probelm: if the window reach the end and h is taking over. It could
% happen that e.g. the h=1 value is empty, that is skiped, but than with
% h=2 the index is subtracted 2 (h=2) which will be empty. The new value is
% therfore calculated from win and plus one empty value.
% Is that realy a problem? It could also just be fine. The value is then
% calculated by less values, still kind of correct.

  
        
%%%%%%%%%% removing nans. Due to the RR distance calculation the first value in NAN
% if exist('RR','var') == 1
%    RR(any(isnan(RR)))=[]; %removing nans
% end
% 

%%%%%%%%%%% CREATING ECG (5min) WINDOWS
 % 300s*500Hz is 5min of data
win_jumps=factor*500;
Fenster=win*500;


if exist('ECG','var') == 1
    m=1;
    uebrig=length(ECG);  % how many minutes are left           

   for k=1:win_jumps:length(ECG)
       if k+Fenster<length(ECG) 
        ECG_win_300{1,m}=ECG(k:k+Fenster-1,1); 
        t_300{1,m}=t_ECG(k:k+Fenster-1,1);
       elseif k+Fenster>=length(ECG) && win_jumps<=uebrig && k>win_jumps*(Fenster-(uebrig/win_jumps)*win_jumps)/win_jumps
%            uebrig=(k+Fenster)-length(ECG);  % how many minutes are left           
           rechts=uebrig/win_jumps;% How many epochs are still left 
           links=(Fenster-rechts*win_jumps)/win_jumps; % how many epochs do we have to atache from the left to get a full 300s window
           ECG_win_300{1,m}=ECG(k-win_jumps*links:k+win_jumps*rechts,1);
           t_300{1,m}=t_ECG(k-win_jumps*links:k+win_jumps*rechts,1);
       elseif k+Fenster>=length(ECG) && win_jumps>uebrig&&  k>win_jumps*(Fenster-(uebrig/win_jumps)*win_jumps)/win_jumps
           rechts=uebrig/win_jumps;% How many epochs are still left 
           links=(Fenster-rechts*win_jumps)/win_jumps; % how many epochs do we have to atache from the left to get a full 300s window
           ECG_win_300{1,m}=ECG(k-win_jumps*links:k+win_jumps*rechts,1);
           t_300{1,m}=t_ECG(k-win_jumps*links:end,1);           
       else
           ECG_win_300{1,m}=ECG(k:end,1);
           t_300{1,m}=t_ECG(k:end,1);           
           
%            break       % if you want to end with the same length for the last cell elementas the others use break. But than the ECG_win_300 is one element shorter thatn ECG_win_30    
       end
       uebrig=length(ECG)-(k+win_jumps);  % how many minutes are left           
       m=m+1;
   end
end

%original&&&&&&&&&&&&&&&&&&&&&&&&& 
%WOrks but does not have same cel size at
%the end
% if exist('ECG','var') == 1
%     m=1;
%    for k=1:win_jumps:length(ECG)
%        if k+Fenster<length(ECG) 
%         ECG_win_300{1,m}=ECG(k:k+Fenster-1,1); 
%         t_300{1,m}=t_ECG(k:k+Fenster-1,1);
%        elseif k+Fenster>=length(ECG) 
%            ECG_win_300{1,m}=ECG(k:end,1);
%            t_300{1,m}=t_ECG(k:end,1);
%            break
%        end
%        m=m+1;
%    end
% end
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


%%%%%%%%%%% CREATING ECG (30s) WINDOWS
 % 300s*500Hz is 5min of data
win_jumps=factor*500;
Fenster=30*500;

if exist('ECG','var') == 1
    m=1;
   for k=1:win_jumps:length(ECG)
       if k+Fenster<length(ECG) 
        ECG_win_30{1,m}=ECG(k:k+Fenster-1,1);
        t_30{1,m}=t_ECG(k:k+Fenster-1,1);        
       elseif k+Fenster>=length(ECG) 
           ECG_win_30{1,m}=ECG(k:end,1);
           t_30{1,m}=t_ECG(k:end,1);           
           break
       end
       m=m+1;
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% SAVING

%  folder=('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\ECG_30secepochs\ECG clustered by annotations\unsorted\');

%%%%%%%%%%%% savig AS QS
% if saving
%    if exist('ECG_win','var') ==1
%        save([folder 'ECG_' num2str(Neonate) '_win_' num2str((factor+1)*30) '.mat'],'ECG_win')  ; 
%    end              
% end % if saving
if saving
    Saving(ECG_win_30,savefolder, Neonate,S)
end
end%win 
%% Nested saving
    function Saving(Feature,savefolder, Neonate, win,S)
        if exist('Feature','var')==1
            name=inputname(1); % variable name of function input
            save([savefolder name '_Session_' num2str(S) '_' num2str(Neonate)],'Feature')
        else
            disp(['saving of ' name ' not possible'])
        end       
    end
