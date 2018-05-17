function SDDec(RR,Neonate,saving,savefolder,win,S) 
%Input
% RR: 5min RR distance data
% Neonate: Which patient
% saving: If saving is whished
% savefolder: Where to save
% win: Duration of the HRV window. Common is 5min/300s

for i=1:length(RR)
    if all(isnan(RR{i,1}))
        SDDEC{1,i}=nan;
        continue
    end    
    [~, columns]=find(RR{i,1}>nanmean(RR{i,1}));
    SDDEC{1,i}=nanstd(RR{i,1}(1,columns));
end
    
%%%%%%%%%%%%replace [] with nan
ix=cellfun(@isempty,SDDEC);
SDDEC(ix)={0};  

            
%%%%%%%%%%%% SAVING            
if saving                     %saving R peaks positions in mat file                 
    Saving(SDDEC,savefolder,Neonate,win,S) 
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
 
 

