function linelength(ECG,t_ECG,Neonate,saving,savefolder,win,Session,S) 
%Input
% RR: 5min RR distance data
% Neonate: Which patient
% saving: If saving is whished
% savefolder: Where to save
% win: Duration of the HRV window. Comon is 5min/300s

[b,a] = butter(3,[0.09 0.29],'stop'); %bandstop with [] in Normalized Frequency x Pi rad/sample
for L=1:length(ECG)
    ECG{1,L} = filtfilt(b,a,ECG{1,L}); %zero phase filter to get rid of QRS complex
    [lineLength{1,L}] = arclength(t_ECG{1,L},(ECG{1,L}),'linear');
end
lineLength{1,1}=lineLength{1,2}; % due to nan, cannot be determined correctly. leaving nan ut still gives false result. rather copy next as they might be similar

%%%%%%%%%%%% SAVING            
    if saving                     %saving R peaks positions in mat file                 
       Saving(lineLength,savefolder,Neonate,win,Session,S) 
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
