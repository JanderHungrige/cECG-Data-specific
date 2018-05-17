function freqdomainHRV (powerspectrum,f,Neonate,win,saving,savefolder,Session,S)



% frequency domain measures of HRV
% this functioncalculates the frequenc donain features over 30 sec epochs,
% while " freqdomain_single.m"  alculates the features over a non splittet 
% signal,ecept for the states AS,QS,...

%%%%%short term analysis (5min)
% 5min total power (<0.4 Hz) % new (<1.5 Hz)
% VLF power in very low frequency range (adult 0.003-0.04 Hz)
% LF  power in low frequency range (adult 0.04-0.15 Hz)
% LFnorm  LF power in normalized units LF/(total power-vlf)*100 or LF/(LF+HF)*100
% HF power in high frequency range (adult 0.15-0.4 Hz)
% HFnorm HF power in normalized units HF/(total power-vlf)*100 or HF/(LF+HF)*100
% LF/HF Ratio LF/HF

%%%%long term analysis (best 24h) [look at function freqdomainHRV_single]
% total power variance of all NN intervals  (<0.4 Hz)
% ULF power in ultra low frequency range (adults <0.003 Hz)
% VLF power in very low frequency range (adult 0.003-0.04 Hz)
% LF  power in low frequency range (adult 0.04-0.15 Hz)
% HF power in high frequency range (adult 0.15-0.4 Hz)
% alpha SLope linear interpolation of the spectrum in a log-log scale (<0.04 Hz)


 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    %% total power (< 1.5 Hz)

upperboundary=1.5;
totpow=zeros(1,length(powerspectrum)); %preallocation
for l=1:length(powerspectrum)
    if isempty(powerspectrum{l,1})==0
        totpow(1,l)=sum(powerspectrum{l,1}(f{1,l}<upperboundary));
        totpow(totpow==0)=nan; %all zeroes to nan to avoid confusion between AS and QS
    end
end
if saving
    Saving(totpow,savefolder, Neonate, win,Session,S)
end
disp(' -totpow calculated')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% VLF power in very low frequency range (adult 0.003-0.04 Hz)

upperboundary=0.04;
lowerboundary=0.003;

VLF=cell(1,length(powerspectrum)); %preallocation
for l=1:length(powerspectrum)
    if isempty(powerspectrum{l,1})==0
        F=f{1,l}<upperboundary;
        F(f{1,l}<lowerboundary)=0;
        VLF_band{1,l}=powerspectrum{l,1}(F,1);%power
    end
end
VLF=cellfun(@sum,VLF_band,'UniformOutput',false);
VLF(find([VLF{:}] == 0))={nan};  % 0 into nan
if saving
    Saving(VLF,savefolder, Neonate, win,Session,S)
end
clearvars F
disp(' -VLF calculated')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
    %% LF  power in low frequency range (adult 0.04-0.15 Hz)
 
upperboundary=0.15;
lowerboundary=0.04;

LF=cell(1,length(powerspectrum)); %preallocation
for l=1:length(powerspectrum)
    if isempty(powerspectrum{l,1})==0
        F=f{1,l}<upperboundary;
        F(f{1,l}<lowerboundary)=0;
        LF_band{1,l}=powerspectrum{l,1}(F,1);%power
    end
end
LF=cellfun(@sum,LF_band,'UniformOutput',false);
LF(find([LF{:}] == 0))={nan};  % 0 into nan
if saving
    Saving(LF,savefolder, Neonate, win,Session,S)
end
clearvars F
disp(' -LF calculated')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% LFnorm  LF power in normalized units LF/(total power-vlf)*100 or LF/(LF+HF)*100


for i=1:length(LF)
    LFnorm(i)=((cell2mat(LF(1,i)))./(totpow(1,i)-cell2mat(VLF(1,i))))*100; % before LFnorm_AS=((cell2mat(LF_AS))./totpowAS)*100;
end
if saving
    Saving(LFnorm,savefolder, Neonate, win,Session,S)
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HF power in high frequency range (adult 0.15-0.4 Hz)

upperboundary=0.4;
lowerboundary=0.15;
HF=cell(1,length(powerspectrum)); %preallocation
for l=1:length(powerspectrum)
    if isempty(powerspectrum{l,1})==0
        F=f{1,l}<upperboundary;
        F(f{1,l}<lowerboundary)=0;
        HF_band{1,l}=powerspectrum{l,1}(F,1);% power
    end
end
HF=cellfun(@sum,HF_band,'UniformOutput',false);
HF(find([HF{:}] == 0))={nan};  % 0 into nan
if saving
    Saving(HF,savefolder, Neonate, win,Session,S)
end
clearvars F
disp(' -HF calculated')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HFnorm HF power in normalized units HF/(total power-vlf)*100 or HF/(LF+HF)*100
    
for i=1:length(HF)
    HFnorm(i)=((cell2mat(HF(1,i)))./(totpow(1,i)-cell2mat(VLF(1,i))))*100;     
end
if saving
    Saving(HFnorm,savefolder, Neonate, win,Session,S)
end
disp(' -HFnorm calculated')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% LF/HF Ratio LF/HF 
    

ratioLFHF=cellfun(@(LF,HF) (LF)/(HF), LF,HF);
if saving
    Saving(ratioLFHF,savefolder, Neonate, win,Session,S)
end   
disp(' -ratioLFHF calculated')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  sHF power in super high frequency range (preterm 0.4-0.7 Hz)

upperboundary=0.7;
lowerboundary=0.4;
sHF=cell(1,length(powerspectrum)); %preallocation
for l=1:length(powerspectrum)
    if isempty(powerspectrum{l,1})==0
        sF=f{1,l}<upperboundary;
        sF(f{1,l}<lowerboundary)=0;
        sHF_band{1,l}=powerspectrum{l,1}(sF,1);% power
    end
end
sHF=cellfun(@sum,sHF_band,'UniformOutput',false);
sHF(find([sHF{:}] == 0))={nan};  % 0 into nan
if saving
    Saving(sHF,savefolder, Neonate, win,Session,S)
end
clearvars F
disp(' -sHF calculated')
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% sHFnorm sHF power in normalized units HF/(total power-vlf)*100 or HF/(LF+HF)*100
    
for i=1:length(sHF)
    sHFnorm(i)=((cell2mat(sHF(1,i)))./(totpow(1,i)-cell2mat(VLF(1,i))))*100;     
end
if saving
    Saving(sHFnorm,savefolder, Neonate, win,Session,S)
end
disp(' -shFnorm calculated')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  uHF power in ultra high frequency range (preterm 0.7-1.5 Hz)

upperboundary=1.5;
lowerboundary=0.7;
uHF=cell(1,length(powerspectrum)); %preallocation
for l=1:length(powerspectrum)
    if isempty(powerspectrum{l,1})==0
        uF=f{1,l}<upperboundary;
        uF(f{1,l}<lowerboundary)=0;
        uHF_band{1,l}=powerspectrum{l,1}(uF,1);% power
    end
end
uHF=cellfun(@sum,uHF_band,'UniformOutput',false);
uHF(find([uHF{:}] == 0))={nan};  % 0 into nan
if saving
    Saving(uHF,savefolder, Neonate, win,Session,S)
end
disp(' -uHF calculated')


clearvars F

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% uHFnorm uHF power in normalized units uHF/(total power-vlf)*100 or uHF/(LF+HF+sHF+uHF)*100
    
for i=1:length(uHF)
    uHFnorm(i)=((cell2mat(uHF(1,i)))./(totpow(1,i)-cell2mat(VLF(1,i))))*100;     
end
if saving
    Saving(uHFnorm,savefolder, Neonate, win,Session,S)
end
disp(' -uHFnorm calculated')



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
 