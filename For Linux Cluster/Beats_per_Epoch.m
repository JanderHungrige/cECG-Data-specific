function Beats_per_Epoch(RR,Neonate,saving,savefolder,win,Session,S)
%Input
% RR: 5min RR distance data
% Neonate: Which patient
% saving: If saving is whished
% savefolder: Where to save
% win: Duration of the HRV window. Comon is 5min/300s

 

    
%%%%%%%%%%%% Calc BpE    

%     BpE=cell(3,length(RR));
    for K=1:length(RR)
        if isempty(RR{K,1})==0 & all(isnan(RR{K,1}))==0
            BpE{1,K}=length(RR{K,1})-1; %-1 as first is nan
        else
            BpE{1,K}=NaN;
        end
    end
   
%%%%%%%%%%%% SAVING            
    if saving                     %saving R peaks positions in mat file                 
       Saving(BpE,savefolder,Neonate,win,Session,S) 
    end% end if saving  
    
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