%The differnce her is that HFnorm and LF norm are not divided by (totpow-VLF) but only by totpow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function freqdomainHRV_fake (pat,plotting,win,saving)
Neonate=pat;

        for j=1:length(win)
                
                
%%%%%%%%%%%%%%%%%%%%%%%%%%short term analysis (5min)
                
%Loading power spectrum and frequency for AS and QS
                        
                cd('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\ECG_30secepochs\spectrum\')

                if exist(fullfile(cd, ['Powerspectrum30_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']), 'file') == 2 % ==> 0 or 2
                    load(fullfile(cd, ['Powerspectrum30_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']))
                end
                if exist(fullfile(cd, ['PowerspectrumAalertness_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']), 'file') == 2 % ==> 0 or 2
                    load(fullfile(cd, ['PowerspectrumAalertness_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']))
                end    
                if exist(fullfile(cd, ['PowerspectrumQalertness_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']), 'file') == 2 % ==> 0 or 2
                    load(fullfile(cd, ['PowerspectrumQalertness_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']))
                end             
                if exist(fullfile(cd, ['PowerspectrumTransition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']), 'file') == 2 % ==> 0 or 2
                    load(fullfile(cd, ['PowerspectrumTransition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']))
                end   
                if exist(fullfile(cd, ['PowerspectrumPosition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']), 'file') == 2 % ==> 0 or 2
                    load(fullfile(cd, ['PowerspectrumPosition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']))
                end    
                if exist(fullfile(cd, ['PowerspectrumNot_reliable_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']), 'file') == 2 % ==> 0 or 2
                    load(fullfile(cd, ['PowerspectrumNot_reliable_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat']))
                end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    %% total power (<0.4Hz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%% QS AS
            if exist('pxxAS','var')
                totpowAS=zeros(1,length(pxxAS)); %preallocation
                for l=1:length(pxxAS)
                    if isempty(pxxAS{1,l})==0
                        totpowAS(1,l)=sum(pxxAS{1,l}(fAS{1,l}<0.4));
                        totpowAS(totpowAS==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
                    end
                end
            end
            
            if exist('pxxQS','var')
                totpowQS=zeros(1,length(pxxQS)); %preallocation                
                for l=1:length(pxxQS)
                    if isempty(pxxQS{1,l})==0
                        totpowQS(1,l)=sum(pxxQS{1,l}(fQS{1,l}<0.4));
                        totpowQS(totpowQS==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
                    end
                end
            end
%%%%%%%%%%%%%% WAKE            
            if exist('pxxAalertness','var')
                totpowAalertness=zeros(1,length(pxxAalertness)); %preallocation
                for l=1:length(pxxAalertness)
                    if isempty(pxxAalertness{1,l})==0
                        totpowAalertness(1,l)=sum(pxxAalertness{1,l}(fAalertness{1,l}<0.4));
                        totpowAalertness(totpowAalertness==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
                    end
                end
            end
            
            if exist('pxxQalertness','var')
                totpowQalertness=zeros(1,length(pxxQalertness)); %preallocation                
                for l=1:length(pxxQalertness)
                    if isempty(pxxQalertness{1,l})==0
                        totpowQalertness(1,l)=sum(pxxQalertness{1,l}(fQalertness{1,l}<0.4));
                        totpowQalertness(totpowQalertness==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
                    end
                end
            end  
%%%%%%%%%%%%%% TRANSITION           
            if exist('pxxtransition','var')
                totpowtransition=zeros(1,length(pxxtransition)); %preallocation
                for l=1:length(pxxtransition)
                    if isempty(pxxtransition{1,l})==0
                        totpowtransition(1,l)=sum(pxxtransition{1,l}(ftransition{1,l}<0.4));
                        totpowtransition(totpowtransition==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
                    end
                end
            end            
%%%%%%%%%%%%%% POSITION          
            if exist('pxxposition','var')
                totpowposition=zeros(1,length(pxxposition)); %preallocation
                for l=1:length(pxxposition)
                    if isempty(pxxposition{1,l})==0
                        totpowposition(1,l)=sum(pxxposition{1,l}(fposition{1,l}<0.4));
                        totpowposition(totpowposition==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
                    end
                end
            end   
%%%%%%%%%%%%%% NOT RELIABLE          
            if exist('pxxNot_reliable','var')
                totpowNot_reliable=zeros(1,length(pxxNot_reliable)); %preallocation
                for l=1:length(pxxNot_reliable)
                    if isempty(pxxNot_reliable{1,l})==0
                        totpowNot_reliable(1,l)=sum(pxxNot_reliable{1,l}(fNot_reliable{1,l}<0.4));
                        totpowNot_reliable(totpowNot_reliable==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
                    end
                end
            end     
%%%%%%%%%%%%%% SAVING            
            if saving
%%%%%% SAVING AS QS                
                folder='C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\HRV features\';
                if exist('totpowAS','var')==1 & exist ('totpowQS','var')
                 save([folder 'totpower_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'totpowAS','totpowQS');
                elseif exist('totpowAS','var') ==1
                 save([folder 'totpower_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'totpowAS');
                elseif exist('totpowQS','var') ==1
                 save([folder 'totpower_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'totpowQS');
                end
%%%%%%% SAVING WAKE                
                if exist('totpowAalertness','var')==1 & exist ('totpowQalertness','var')
                 save([folder 'totpowerWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'totpowAalertness','totpowQalertness');
                elseif exist('totpowAalertness','var') ==1
                 save([folder 'totpowerWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'totpowAalertness');
                elseif exist('totpowQalertness','var') ==1
                 save([folder 'totpowerWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'totpowQalertness');
                end
%%%%%% SAVING TRANSITION
                if exist('totpowtransition','var')==1
                 save([folder 'totpowertransition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'totpowtransition');
                end
%%%%%% SAVING POSITION
                if exist('totpowposition','var')==1
                 save([folder 'totpowposition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'totpowposition');
                end 
%%%%%% SAVING POSITION
                if exist('totpowNot_reliable','var')==1
                 save([folder 'totpowNot_reliable_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'totpowNot_reliable');
                end                 
             
            end% if saving
Berechnung7=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% VLF power in very low frequency range (adult 0.003-0.04 Hz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% AS QS
            if exist('pxxAS','var')
                VLF_AS=cell(1,length(pxxAS)); %preallocation
                for l=1:length(pxxAS)
                    if isempty(pxxAS{1,l})==0
                        FAS=fAS{1,l}<0.04;
                        FAS(fAS{1,l}<0.003)=0;
%                         VLF_AS_band{1,l}=fAS{1,l}(FAS,1); % determine the frequency bands
                        VLF_AS_band{1,l}=pxxAS{1,l}(FAS,1);%power
                    end
                end
                VLF_AS=cellfun(@sum,VLF_AS_band,'UniformOutput',false);
                VLF_AS(find([VLF_AS{:}] == 0))={nan};  % 0 into nan
            end

            if exist('pxxQS','var')
                VLF_QS=cell(1,length(pxxQS)); %preallocation
                for l=1:length(pxxQS)
                    if isempty(pxxQS{1,l})==0
                        FQS=fQS{1,l}<0.04;
                        FQS(fQS{1,l}<0.003)=0;
%                         VLF_QS_band{1,l}=fQS{1,l}(FQS);%frequency
                        VLF_QS_band{1,l}=pxxQS{1,l}(FQS); %power
                    end
                end
                VLF_QS=cellfun(@sum,VLF_QS_band,'UniformOutput',false);
                VLF_QS(find([VLF_QS{:}] == 0))={nan};  % 0 into nan
            end
%%%%%%%% WAKE            
            if exist('pxxAalertness','var')
                VLF_Aalertness=cell(1,length(pxxAalertness)); %preallocation
                for l=1:length(pxxAalertness)
                    if isempty(pxxAalertness{1,l})==0
                        FAalertness=fAalertness{1,l}<0.04;
                        FAalertness(fAalertness{1,l}<0.003)=0;
%                         VLF_AS_band{1,l}=fAS{1,l}(FAS,1); % determine the frequency bands
                        VLF_Aalertness_band{1,l}=pxxAalertness{1,l}(FAalertness,1);%power
                    end
                end
                VLF_Aalertness=cellfun(@sum,VLF_Aalertness_band,'UniformOutput',false);
                VLF_Aalertness(find([VLF_Aalertness{:}] == 0))={nan};  % 0 into nan
            end

            if exist('pxxQalertness','var')
                VLF_Qalertness=cell(1,length(pxxQalertness)); %preallocation
                for l=1:length(pxxQalertness)
                    if isempty(pxxQalertness{1,l})==0
                        FQalertness=fQalertness{1,l}<0.04;
                        FQalertness(fQalertness{1,l}<0.003)=0;
%                         VLF_QS_band{1,l}=fQS{1,l}(FQS);%frequency
                        VLF_Qalertness_band{1,l}=pxxQalertness{1,l}(FQalertness); %power
                    end
                end
                VLF_Qalertness=cellfun(@sum,VLF_Qalertness_band,'UniformOutput',false);
                VLF_Qalertness(find([VLF_Qalertness{:}] == 0))={nan};  % 0 into nan
            end  
%%%%%%%% TRANSITION       
            if exist('pxxtransition','var')
                VLF_transition=cell(1,length(pxxtransition)); %preallocation
                for l=1:length(pxxtransition)
                    if isempty(pxxtransition{1,l})==0
                        Ftransition=ftransition{1,l}<0.04;
                        Ftransition(ftransition{1,l}<0.003)=0;
%                         VLF_AS_band{1,l}=fAS{1,l}(FAS,1); % determine the frequency bands
                        VLF_transition_band{1,l}=pxxtransition{1,l}(Ftransition,1);%power
                    end
                end
                VLF_transition=cellfun(@sum,VLF_transition_band,'UniformOutput',false);
                VLF_transition(find([VLF_transition{:}] == 0))={nan};  % 0 into nan
            end           
%%%%%%%% POSITION       
            if exist('pxxposition','var')
                VLF_position=cell(1,length(pxxposition)); %preallocation
                for l=1:length(pxxposition)
                    if isempty(pxxposition{1,l})==0
                        Fposition=fposition{1,l}<0.04;
                        Fposition(fposition{1,l}<0.003)=0;
%                         VLF_AS_band{1,l}=fAS{1,l}(FAS,1); % determine the frequency bands
                        VLF_position_band{1,l}=pxxposition{1,l}(Fposition,1);%power
                    end
                end
                VLF_position=cellfun(@sum,VLF_position_band,'UniformOutput',false);
                VLF_position(find([VLF_position{:}] == 0))={nan};  % 0 into nan
            end 
%%%%%%%% NOT RELIABLE       
            if exist('pxxNot_reliable','var')
                VLF_Not_reliable=cell(1,length(pxxNot_reliable)); %preallocation
                for l=1:length(pxxNot_reliable)
                    if isempty(pxxNot_reliable{1,l})==0
                        FNot_reliable=fNot_reliable{1,l}<0.04;
                        FNot_reliable(fNot_reliable{1,l}<0.003)=0;
%                         VLF_AS_band{1,l}=fAS{1,l}(FAS,1); % determine the frequency bands
                        VLF_Not_reliable_band{1,l}=pxxNot_reliable{1,l}(FNot_reliable,1);%power
                    end
                end
                VLF_Not_reliable=cellfun(@sum,VLF_Not_reliable_band,'UniformOutput',false);
                VLF_Not_reliable(find([VLF_Not_reliable{:}] == 0))={nan};  % 0 into nan
            end             
            
%%%%%%%%%% SAVING            
            if saving
                folder='C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\HRV features\';
                
%%%%% SAVING AS QS                
                if exist('VLF_AS','var')==1 & exist ('VLF_QS','var')
                 save([folder 'VLF_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'VLF_AS','VLF_QS');
                elseif exist('VLF_AS','var') ==1
                 save([folder 'VLF_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'VLF_AS');
                elseif exist('VLF_QS','var') ==1
                 save([folder 'VLF_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'VLF_QS');
                end
%%%%% SAVING WAKE                
                if exist('VLF_Aalertness','var')==1 & exist ('VLF_Qalertness','var')
                 save([folder 'VLFWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'VLF_Aalertness','VLF_Qalertness');
                elseif exist('VLF_Aalertness','var') ==1
                 save([folder 'VLFWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'VLF_Aalertness');
                elseif exist('VLF_Qalertness','var') ==1
                 save([folder 'VLFWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'VLF_Qalertness');
                end  
%%%%% SAVING TRANSITION              
                if exist('VLF_transition','var')==1
                 save([folder 'VLFtransition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'VLF_transition');                
                end
%%%%% SAVING POSITION              
                if exist('VLF_position','var')==1
                 save([folder 'VLFposition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'VLF_position');                
                end      
%%%%% SAVING NOT RELIABLE              
                if exist('VLF_Not_reliable','var')==1
                 save([folder 'VLFNot_reliable_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'VLF_Not_reliable');                
                end                  
            
            end % IF SAVING

                    clearvars FQS FAS
Berechnung7=2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
    %% LF  power in low frequency range (adult 0.04-0.15 Hz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%% AS QS
            if exist('pxxAS','var')
                LF_AS=cell(1,length(pxxAS)); %preallocation
                for l=1:length(pxxAS)
                    if isempty(pxxAS{1,l})==0
                        FAS=fAS{1,l}<0.15;
                        FAS(fAS{1,l}<0.04)=0;
%                         LF_AS_band{1,l}=fAS{1,l}(FAS,1);%frequency
                        LF_AS_band{1,l}=pxxAS{1,l}(FAS,1);%power
                    end
                end
                LF_AS=cellfun(@sum,LF_AS_band,'UniformOutput',false);
                LF_AS(find([LF_AS{:}] == 0))={nan};  % 0 into nan
            end

            if exist('pxxQS','var')
                LF_QS=cell(1,length(pxxQS)); %preallocation
                for l=1:length(pxxQS)
                    if isempty(pxxQS{1,l})==0
                        FQS=fQS{1,l}<0.15;
                        FQS(fQS{1,l}<0.04)=0;
%                         LF_QS_band{1,l}=fQS{1,l}(FQS);% frequency 
                        LF_QS_band{1,l}=pxxQS{1,l}(FQS);% power 
                    end
                end
                LF_QS=cellfun(@sum,LF_QS_band,'UniformOutput',false);
                LF_QS(find([LF_QS{:}] == 0))={nan};  % 0 into nan
            end
%%%%% WAKE
            if exist('pxxAalertness','var')
                LF_Aalertness=cell(1,length(pxxAalertness)); %preallocation
                for l=1:length(pxxAalertness)
                    if isempty(pxxAalertness{1,l})==0
                        FAalertness=fAalertness{1,l}<0.15;
                        FAalertness(fAalertness{1,l}<0.04)=0;
%                         LF_AS_band{1,l}=fAS{1,l}(FAS,1);%frequency
                        LF_Aalertness_band{1,l}=pxxAalertness{1,l}(FAalertness,1);%power
                    end
                end
                LF_Aalertness=cellfun(@sum,LF_Aalertness_band,'UniformOutput',false);
                LF_Aalertness(find([LF_Aalertness{:}] == 0))={nan};  % 0 into nan
            end

            if exist('pxxQalertness','var')
                LF_Qalertness=cell(1,length(pxxQalertness)); %preallocation
                for l=1:length(pxxQalertness)
                    if isempty(pxxQalertness{1,l})==0
                        FQalertness=fQalertness{1,l}<0.15;
                        FQalertness(fQalertness{1,l}<0.04)=0;
%                         LF_QS_band{1,l}=fQS{1,l}(FQS);% frequency 
                        LF_Qalertness_band{1,l}=pxxQalertness{1,l}(FQalertness);% power 
                    end
                end
                LF_Qalertness=cellfun(@sum,LF_Qalertness_band,'UniformOutput',false);
                LF_Qalertness(find([LF_Qalertness{:}] == 0))={nan};  % 0 into nan
            end  
%%%%% TRANSITION
            if exist('pxxtransition','var')
                LF_transition=cell(1,length(pxxtransition)); %preallocation
                for l=1:length(pxxtransition)
                    if isempty(pxxtransition{1,l})==0
                        Ftransition=ftransition{1,l}<0.15;
                        Ftransition(ftransition{1,l}<0.04)=0;
%                         LF_AS_band{1,l}=fAS{1,l}(FAS,1);%frequency
                        LF_transition_band{1,l}=pxxtransition{1,l}(Ftransition,1);%power
                    end
                end
                LF_transition=cellfun(@sum,LF_transition_band,'UniformOutput',false);
                LF_transition(find([LF_transition{:}] == 0))={nan};  % 0 into nan
            end  
%%%%% POSITION
            if exist('pxxposition','var')
                LF_position=cell(1,length(pxxposition)); %preallocation
                for l=1:length(pxxposition)
                    if isempty(pxxposition{1,l})==0
                        Fposition=fposition{1,l}<0.15;
                        Fposition(fposition{1,l}<0.04)=0;
%                         LF_AS_band{1,l}=fAS{1,l}(FAS,1);%frequency
                        LF_position_band{1,l}=pxxposition{1,l}(Fposition,1);%power
                    end
                end
                LF_position=cellfun(@sum,LF_position_band,'UniformOutput',false);
                LF_position(find([LF_position{:}] == 0))={nan};  % 0 into nan
            end     
%%%%% NOT RELIABLE
            if exist('pxxNot_reliable','var')
                LF_Not_reliable=cell(1,length(pxxNot_reliable)); %preallocation
                for l=1:length(pxxNot_reliable)
                    if isempty(pxxNot_reliable{1,l})==0
                        FNot_reliable=fNot_reliable{1,l}<0.15;
                        FNot_reliable(fNot_reliable{1,l}<0.04)=0;
%                         LF_AS_band{1,l}=fAS{1,l}(FAS,1);%frequency
                        LF_Not_reliable_band{1,l}=pxxNot_reliable{1,l}(FNot_reliable,1);%power
                    end
                end
                LF_Not_reliable=cellfun(@sum,LF_Not_reliable_band,'UniformOutput',false);
                LF_Not_reliable(find([LF_Not_reliable{:}] == 0))={nan};  % 0 into nan
            end              
%%%%% SAVING
            if saving
                folder='C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\HRV features\';
%%% SAVING AS QS                
                if exist('LF_AS','var' )==1 & exist ('LF_QS','var')
                 save([folder 'LF_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'LF_AS','LF_QS');
                elseif exist('LF_AS','var') ==1
                 save([folder 'LF_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'LF_AS');
                elseif exist('LF_QS','var') ==1
                 save([folder 'LF_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'LF_QS');
                end

%%% SAVING WAKE
                if exist('LF_Aalertness','var')==1 & exist ('LF_Qalertness','var')
                 save([folder 'LFWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'LF_Aalertness','LF_Qalertness');
                elseif exist('LF_Aalertness','var') ==1
                 save([folder 'LFWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'LF_Aalertness');
                elseif exist('LF_Qalertness','var') ==1
                 save([folder 'LFWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'LF_Qalertness');
                end
%%% SAVING TRANSITION
                if exist('LF_transition','var')==1 
                 save([folder 'LFtransition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'LF_transition');                
                end
%%% SAVING POSITION
                if exist('LF_position','var')==1 
                 save([folder 'LFposition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'LF_position');                
                end           
%%% SAVING NOT RELIABLE
                if exist('LF_Not_reliable','var')==1 
                 save([folder 'LF_Not_reliable_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'LF_Not_reliable');                
                end                   
             
            
            end% if saving
                    clearvars FQS FAS
Berechnung7=3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% LFnorm  LF power in normalized units LF/(total power-vlf)*100 or LF/(LF+HF)*100
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% AS QS
          if exist('pxxAS','var')
                    LFnorm_AS=((cell2mat(LF_AS))./(totpowAS))*100; % before LFnorm_AS=((cell2mat(LF_AS))./totpowAS)*100;
                    LFnorm_AS(LFnorm_AS==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
          end
          
          if exist('pxxQS','var')
                    LFnorm_QS=((cell2mat(LF_QS))./(totpowQS))*100;
                    LFnorm_QS(LFnorm_QS==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
          end
%%%%% WAKE
          if exist('pxxAalertness','var')
                    LFnorm_Aalertness=((cell2mat(LF_Aalertness))./(totpowAalertness-cell2mat(VLF_Aalertness)))*100; 
                    LFnorm_Aalertness(LFnorm_Aalertness==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
          end
          
          if exist('pxxQalertness','var')
                    LFnorm_Qalertness=((cell2mat(LF_Qalertness))./(totpowQalertness-cell2mat(VLF_Qalertness)))*100;
                    LFnorm_Qalertness(LFnorm_Qalertness==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
          end
%%%%% TRANSITION
          if exist('pxxtransition','var')
                    LFnorm_transition=((cell2mat(LF_transition))./(totpowtransition-cell2mat(VLF_transition)))*100; 
                    LFnorm_transition(LFnorm_transition==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
          end     
%%%%% POSITION
          if exist('pxxposition','var')
                    LFnorm_position=((cell2mat(LF_position))./(totpowposition-cell2mat(VLF_position)))*100; 
                    LFnorm_position(LFnorm_position==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
          end            
%%%%% NOT RELIABLE
          if exist('pxxNot_reliable','var')
                    LFnorm_Not_reliable=((cell2mat(LF_Not_reliable))./(totpowNot_reliable-cell2mat(VLF_Not_reliable)))*100; 
                    LFnorm_Not_reliable(LFnorm_Not_reliable==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
          end    
%%%% SAVING
            if saving
                folder='C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\HRV features\';
%%%% SAVING AS QS                
                if exist('LFnorm_AS','var')==1 & exist ('LFnorm_QS','var')
                 save([folder 'LFnorm_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'LFnorm_AS','LFnorm_QS');
                elseif exist('LFnorm_AS','var') ==1
                 save([folder 'LFnorm_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'LFnorm_AS');
                elseif exist('LFnorm_QS','var') ==1
                 save([folder 'LFnorm_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'LFnorm_QS');
                end
%%%% SAVING WAKE                
                if exist('LFnorm_Aalertness','var')==1 & exist ('LFnorm_Qalertness','var')
                 save([folder 'LFnormWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'LFnorm_Aalertness','LFnorm_Qalertness');
                elseif exist('LFnorm_Aalertness','var') ==1
                 save([folder 'LFnormWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'LFnorm_Aalertness');
                elseif exist('LFnorm_Qalertness','var') ==1
                 save([folder 'LFnormWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'LFnorm_Qalertness');
                end  
%%%% SAVING TRANSITION                
                if exist('LFnorm_transition','var')==1 
                 save([folder 'LFnormtransition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'LFnorm_transition');                
                end
%%%% SAVING POSITION                
                if exist('LFnorm_position','var')==1 
                 save([folder 'LFnormposition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'LFnorm_position');                
                end  
%%%% SAVING NOT RELIABLE                
                if exist('LFnorm_Not_reliable','var')==1 
                 save([folder 'LFnormNot_reliable_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'LFnorm_Not_reliable');                
                end                                 
            end   % if saving     
Berechnung7=4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HF power in high frequency range (adult 0.15-0.4 Hz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% AS QS
          if exist('pxxAS','var')
                    HF_AS=cell(1,length(pxxAS)); %preallocation
                    for l=1:length(pxxAS)
                        if isempty(pxxAS{1,l})==0
                            FAS=fAS{1,l}<0.4;
                            FAS(fAS{1,l}<0.15)=0;
    %                         HF_AS_band{1,l}=fAS{1,l}(FAS,1);% frequency
                            HF_AS_band{1,l}=pxxAS{1,l}(FAS,1);% power
                        end
                    end
                    HF_AS=cellfun(@sum,HF_AS_band,'UniformOutput',false);
                    HF_AS(find([HF_AS{:}] == 0))={nan};  % 0 into nan
          end

          if exist('pxxQS','var')          
                    HF_QS=cell(1,length(pxxQS)); %preallocation
                    for l=1:length(pxxQS)
                        if isempty(pxxQS{1,l})==0
                            FQS=fQS{1,l}<0.4;
                            FQS(fQS{1,l}<0.15)=0;
%                             HF_QS_band{1,l}=fQS{1,l}(FQS);% frequency
                            HF_QS_band{1,l}=pxxQS{1,l}(FQS);% power
                        end
                    end
                    HF_QS=cellfun(@sum,HF_QS_band,'UniformOutput',false);
                    HF_QS(find([HF_QS{:}] == 0))={nan};  % 0 into nan
          end
%%%%% WAKE
          if exist('pxxAalertness','var')
                    HF_Aalertness=cell(1,length(pxxAalertness)); %preallocation
                    for l=1:length(pxxAalertness)
                        if isempty(pxxAalertness{1,l})==0
                            FAalertness=fAalertness{1,l}<0.4;
                            FAalertness(fAalertness{1,l}<0.15)=0;
    %                         HF_AS_band{1,l}=fAS{1,l}(FAS,1);% frequency
                            HF_Aalertness_band{1,l}=pxxAalertness{1,l}(FAalertness,1);% power
                        end
                    end
                    HF_Aalertness=cellfun(@sum,HF_Aalertness_band,'UniformOutput',false);
                    HF_Aalertness(find([HF_Aalertness{:}] == 0))={nan};  % 0 into nan
          end
          if exist('pxxQalertness','var')          
                    HF_Qalertness=cell(1,length(pxxQalertness)); %preallocation
                    for l=1:length(pxxQalertness)
                        if isempty(pxxQalertness{1,l})==0
                            FQalertness=fQalertness{1,l}<0.4;
                            FQalertness(fQalertness{1,l}<0.15)=0;
%                             HF_Qalertness_band{1,l}=fQalertness{1,l}(FQalertness);% frequency
                            HF_Qalertness_band{1,l}=pxxQalertness{1,l}(FQalertness);% power
                        end
                    end
                    HF_Qalertness=cellfun(@sum,HF_Qalertness_band,'UniformOutput',false);
                    HF_Qalertness(find([HF_Qalertness{:}] == 0))={nan};  % 0 into nan
          end  
%%%%% TRANSITION
          if exist('pxxtransition','var')
                    HF_transition=cell(1,length(pxxtransition)); %preallocation
                    for l=1:length(pxxtransition)
                        if isempty(pxxtransition{1,l})==0
                            Ftransition=ftransition{1,l}<0.4;
                            Ftransition(ftransition{1,l}<0.15)=0;
    %                         HF_AS_band{1,l}=fAS{1,l}(FAS,1);% frequency
                            HF_transition_band{1,l}=pxxtransition{1,l}(Ftransition,1);% power
                        end
                    end
                    HF_transition=cellfun(@sum,HF_transition_band,'UniformOutput',false);
                    HF_transition(find([HF_transition{:}] == 0))={nan};  % 0 into nan
          end   
%%%%% POSITION
          if exist('pxxposition','var')
                    HF_position=cell(1,length(pxxposition)); %preallocation
                    for l=1:length(pxxposition)
                        if isempty(pxxposition{1,l})==0
                            Fposition=fposition{1,l}<0.4;
                            Fposition(fposition{1,l}<0.15)=0;
    %                         HF_AS_band{1,l}=fAS{1,l}(FAS,1);% frequency
                            HF_position_band{1,l}=pxxposition{1,l}(Fposition,1);% power
                        end
                    end
                    HF_position=cellfun(@sum,HF_position_band,'UniformOutput',false);
                    HF_position(find([HF_position{:}] == 0))={nan};  % 0 into nan
          end  
%%%%% NOT RELIABLE
          if exist('pxxNot_reliable','var')
                    HF_Not_reliable=cell(1,length(pxxNot_reliable)); %preallocation
                    for l=1:length(pxxNot_reliable)
                        if isempty(pxxNot_reliable{1,l})==0
                            FNot_reliable=fNot_reliable{1,l}<0.4;
                            FNot_reliable(fNot_reliable{1,l}<0.15)=0;
    %                         HF_Not_reliable_band{1,l}=fNot_reliable{1,l}(FNot_reliable,1);% frequency
                            HF_Not_reliable_band{1,l}=pxxNot_reliable{1,l}(FNot_reliable,1);% power
                        end
                    end
                    HF_Not_reliable=cellfun(@sum,fNot_reliable,'UniformOutput',false);
                    HF_Not_reliable(find([HF_Not_reliable{:}] == 0))={nan};  % 0 into nan
          end           
%%%%%% SAVING
             if saving
                folder='C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\HRV features\';
%%%%% SAVING AS QS                
                if exist('HF_AS','var')==1 & exist ('HF_QS','var')
                 save([folder 'HF_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'HF_AS','HF_QS');
                elseif exist('HF_AS','var') ==1
                 save([folder 'HF_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'HF_AS');
                elseif exist('HF_QS','var') ==1
                 save([folder 'HF_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'HF_QS');
                end             
%%%%% SAVING WAKE                
                if exist('HF_Aalertness','var')==1 & exist ('HF_Qalertness','var')
                 save([folder 'HFWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'HF_Aalertness','HF_Qalertness');
                elseif exist('HF_Aalertness','var') ==1
                 save([folder 'HFWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'HF_Aalertness');
                elseif exist('HF_Qalertness','var') ==1
                 save([folder 'HFWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'HF_Qalertness');
                end          
%%%%% SAVING TRANSITION                
                if exist('HF_transition','var')==1 
                 save([folder 'HFtransition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'HF_transition');             
                end
%%%%% SAVING POSITION                
                if exist('HF_position','var')==1 
                 save([folder 'HFposition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'HF_position');             
                end   
%%%%% SAVING NOT RELIABLE                
                if exist('HF_Not_reliable','var')==1 
                 save([folder 'HFNot_reliable_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'HF_Not_reliable');             
                end   
             end% if saving

                    clearvars FQS FAS
                    Berechnung7=5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HFnorm HF power in normalized units HF/(total power-vlf)*100 or HF/(LF+HF)*100
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%% AS QS
              if exist('pxxAS','var')
                    HFnorm_AS=((cell2mat(HF_AS))./(totpowAS))*100;     
                    HFnorm_AS(HFnorm_AS==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
              end
              
              if exist('pxxQS','var')
                    HFnorm_QS=((cell2mat(HF_QS))./(totpowQS))*100;     
                    HFnorm_QS(HFnorm_QS==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
              end
%%%%%%%%% WAKE
              if exist('pxxAalertness','var')
                    HFnorm_Aalertness=((cell2mat(HF_Aalertness))./(totpowAalertness-cell2mat(VLF_Aalertness)))*100;     
                    HFnorm_Aalertness(HFnorm_Aalertness==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
              end
              
              if exist('pxxQalertness','var')
                    HFnorm_Qalertness=((cell2mat(HF_Qalertness))./(totpowQalertness-cell2mat(VLF_Qalertness)))*100;     
                    HFnorm_Qalertness(HFnorm_Qalertness==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
              end   
%%%%%%%%% TRANSITION
              if exist('pxxtransition','var')
                    HFnorm_transition=((cell2mat(HF_transition))./(totpowtransition-cell2mat(VLF_transition)))*100;     
                    HFnorm_transition(HFnorm_transition==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
              end       
%%%%%%%%% POSITION
              if exist('pxxposition','var')
                    HFnorm_position=((cell2mat(HF_position))./(totpowposition-cell2mat(VLF_position)))*100;     
                    HFnorm_position(HFnorm_position==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
              end                 
%%%%%%%%% NOT RELIABLE
              if exist('pxxNot_reliable','var')
                    HFnorm_Not_reliable=((cell2mat(HF_Not_reliable))./(totpowNot_reliable-cell2mat(VLF_Not_reliable)))*100;     
                    HFnorm_Not_reliable(HFnorm_Not_reliable==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
              end                
%%%%%% SAVING              
             if saving
                folder='C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\HRV features\';
%%%%%% SAVING AS QS                
                if exist('HFnorm_AS','var')==1 & exist ('HFnorm_QS','var')
                 save([folder 'HFnorm_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'HFnorm_AS','HFnorm_QS');
                elseif exist('HFnorm_AS','var') ==1
                 save([folder 'HFnorm_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'HFnorm_AS');
                elseif exist('HFnorm_QS','var') ==1
                 save([folder 'HFnorm_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'HFnorm_QS');
                end
%%%%%% SAVING WAKE                
                if exist('HFnorm_Aalertness','var')==1 & exist ('HFnorm_Qalertness','var')
                 save([folder 'HFnormWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'HFnorm_Aalertness','HFnorm_Qalertness');
                elseif exist('HFnorm_Aalertness','var') ==1
                 save([folder 'HFnormWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'HFnorm_Aalertness');
                elseif exist('HFnorm_Qalertness','var') ==1
                 save([folder 'HFnormWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'HFnorm_Qalertness');
                end 
%%%%%% SAVING TRANSITION                
                if exist('HFnorm_transition','var')==1 
                 save([folder 'HFnormtransition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'HFnorm_transition');
                end
%%%%%% SAVING POSITION                
                if exist('HFnorm_position','var')==1 
                 save([folder 'HFnormposition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'HFnorm_position');
                end   
%%%%%% SAVING NOT RELIABLE                
                if exist('HFnorm_Not_reliable','var')==1 
                 save([folder 'HFnorm_Not_reliable_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'HFnorm_Not_reliable');
                end                  
            
             end % IF SAVING       
Berechnung7=6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% LF/HF Ratio LF/HF 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%% AS QS
              if exist('pxxAS','var')
                    ratioLFHF_AS=cellfun(@(LF_AS,HF_AS) (LF_AS)/(HF_AS), LF_AS,HF_AS);
              end
              if exist('pxxQS','var')              
                    ratioLFHF_QS=cellfun(@(LF_QS,HF_QS) LF_QS/HF_QS, LF_QS,HF_QS);
              end
%%%%% WAKE
              if exist('pxxAalertness','var')
                    ratioLFHF_Aalertness=cellfun(@(LF_Aalertness,HF_Aalertness) LF_Aalertness/HF_Aalertness, LF_Aalertness,HF_Aalertness);
              end
              if exist('pxxQalertness','var')              
                    ratioLFHF_Qalertness=cellfun(@(LF_Qalertness,HF_Qalertness) LF_Qalertness/HF_Qalertness, LF_Qalertness,HF_Qalertness);
              end    
%%%%% TRANSITION
              if exist('pxxtransition','var')
                    ratioLFHF_transition=cellfun(@(LF_transition,HF_transition) LF_transition/HF_transition, LF_transition,HF_transition);
              end  
%%%%% POSITION
              if exist('pxxposition','var')
                    ratioLFHF_position=cellfun(@(LF_position,HF_position) LF_position/HF_position, LF_position,HF_position);
              end
%%%%% NOT RELIABLE
              if exist('pxxNot_reliable','var')
                    ratioLFHF_Not_reliable=cellfun(@(LF_Not_reliable,HF_Not_reliable) LF_Not_reliable/HF_Not_reliable, LF_Not_reliable,HF_Not_reliable);
              end              
%%%%%%% SAVING
             if saving
                folder='C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\HRV features\';
%%%% SAVING AS QS                
                if exist('ratioLFHF_AS','var')==1 & exist ('ratioLFHF_QS','var')
                 save([folder 'ratioLFHF_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'ratioLFHF_AS','ratioLFHF_QS');
                elseif exist('ratioLFHF_AS','var') ==1
                 save([folder 'ratioLFHF_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'ratioLFHF_AS');
                elseif exist('ratioLFHF_QS','var') ==1
                 save([folder 'ratioLFHF_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'ratioLFHF_QS');
                end
%%%% SAVING WAKE               
                if exist('ratioLFHF_Aalertness','var')==1 & exist ('ratioLFHF_Qalertness','var')
                 save([folder 'ratioLFHFWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'ratioLFHF_Aalertness','ratioLFHF_Qalertness');
                elseif exist('ratioLFHF_Aalertness','var') ==1
                 save([folder 'ratioLFHFWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'ratioLFHF_Aalertness');
                elseif exist('ratioLFHF_Qalertness','var') ==1
                 save([folder 'ratioLFHFWake_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'ratioLFHF_Qalertness');
                end       
%%%% SAVING TRANSITION               
                if exist('ratioLFHF_transition','var')==1 
                 save([folder 'ratioLFHFtransition_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'ratioLFHF_transition');                
                end
%%%% SAVING POSITION               
                if exist('ratioLFHF_position','var')==1 
                 save([folder 'ratioLFHF_position_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'ratioLFHF_position');                
                end 
%%%% SAVING NOT RELIABLE               
                if exist('ratioLFHF_Not_reliable','var')==1 
                 save([folder 'ratioLFHF_Not_reliable_' num2str(Neonate) '_win_' num2str(win(1,j)) '.mat'],'ratioLFHF_Not_reliable');                
                end                        
             end   % if saving  

Berechnung7=7
                clearvars -except Neonate win i j k saving plotting Berechnung7
            
        end % win end
end

