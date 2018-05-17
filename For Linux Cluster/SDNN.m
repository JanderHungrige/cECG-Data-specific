function SDNN(RR,Neonate,saving,savefolder,win,Session,S) 
%Input
% RR: 5min RR distance data
% Neonate: Which patient
% saving: If saving is whished
% savefolder: Where to save
% win: Duration of the HRV window. Comon is 5min/300s

for i=1:length(RR)
     SDNN{1,i}=nanstd(RR{i,1});
end


%%%%%%%%%%%%replace [] with nan
ix=cellfun(@isempty,SDNN);
SDNN(ix)={nan};  

            
%%%%%%%%%%%% SAVING            
if saving                     %saving R peaks positions in mat file                 
    Saving(SDNN,savefolder,Neonate,win,Session,S) 
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
 
 
 

