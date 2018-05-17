function ECG_HRV_power(powerspectrum,HRV,ECG,HRV300,ECG300,Neonate,saving,savefolder,win,Session,S) 
%
%ECG/HRV power is 30s win
%%%%%%%%%%%% SAVING     
    for i =1:length(HRV)
        HRV{i,1}=HRV{i,1}' ;
    end
    if saving                     %saving R peaks positions in mat file 
       HRV=HRV';
       Saving(ECG,savefolder,Neonate,30,Session,S) 
       Saving(HRV,savefolder,Neonate,30,Session,S)        
       ECG=ECG300;
       HRV=HRV300' ;
       powerspectrum=powerspectrum' ;
       Saving(ECG,savefolder,Neonate,300,Session,S) 
       Saving(HRV,savefolder,Neonate,300,Session,S)        
       Saving(powerspectrum,savefolder,Neonate,300,Session,S)        
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
