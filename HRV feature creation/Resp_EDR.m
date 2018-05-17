function Resp_EDR(Resp300,Resp30,EDR300,EDR30,Neonate,saving,savefolder,win,Session,S) 
%
%ECG/HRV power is 30s win
%%%%%%%%%%%% SAVING            
    if saving                     %saving R peaks positions in mat file    
       Resp=Resp300;
       EDR=EDR300;
       Saving(Resp,savefolder,Neonate,300,Session,S) 
       Saving(EDR,savefolder,Neonate,300,Session,S) 
       Resp=Resp30;
       EDR=EDR30;
       Saving(EDR,savefolder,Neonate,30,Session,S)   
       Saving(Resp,savefolder,Neonate,30,Session,S) 
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
