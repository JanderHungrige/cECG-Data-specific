function RMSSD(RR,Neonate,saving,savefolder,win,Session,S) 
%Input
% RR: 5min RR distance data
% Neonate: Which patient
% saving: If saving is whished
% savefolder: Where to save
% win: Duration of the HRV window. Comon is 5min/300s

% ouput is the RMSSD, which is the sqrt(MSE) of the 5 minute epochs of the NN
% intervales. 

%deleting NaNs
idx=cellfun(@isnan,RR,'UniformOutput',false);
for i=1:length(RR)
    RR{i,1}(idx{i,1})=[]; 
end


for i=1:length(RR)
    if isempty(RR{i,1})==0
        RR{i,1}(2,:)=circshift(RR{i,1}',[1 0]);   %shifting the array to gain a Y and a Y+1(estimated) with wich we can determine the MSE
        [PSNR,MSE{1,i},MAXERR,L2RAT] = measerr(RR{i,1}(1),RR{i,1}(2)); %Calculating MSE
        RMSSD=cellfun(@sqrt,MSE, 'Un', false); %sqare root over the MSE of the adjacent NN intervals creates the RMSSD
    end
end

%%%%%%%%%%%%replace [] with nan   

ix=cellfun(@isempty,RMSSD);
RMSSD(ix)={nan};  


%%%%%%%%%%%% SAVING            
if saving                     %saving R peaks positions in mat file                 
   Saving(RMSSD,savefolder,Neonate,win,Session,S) 
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
 
 