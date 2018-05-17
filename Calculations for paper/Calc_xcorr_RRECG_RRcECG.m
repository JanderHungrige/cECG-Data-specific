
%Load ECG and cECG data ********************************************************
addpath('C:\Users\310122653\Documents\PhD\Matlab\R peak detection and HRV')
addpath('C:\Users\310122653\Documents\PhD\Matlab\cECG Data specific\R peak detection')
datafolder='E:\cECG_study\B_Annotations\Datafiles\For Quick Annotator';
Methode='M';   % R for Ralph; M for Michiel
Sessiondir=dir([datafolder '\participant*' ]);

Patients=[5,6,7,9,10,11,12];

Ralph=[1;-1;1; 1; 1;-1; 1];%Determine if the ECG signal should be turned -1 or not 1. 
      %4 5  6 7  9  10 11 12 13
Corr_S_Cell=cell(length(Patients),1);
Corr_Session=[];
Corr_Patient=[];

for P=1:length(Patients)
    Sessiondir=dir([datafolder '\participant' num2str(Patients(P)) '\DAQ\2012*' ]);
    Sessions=dir([Sessiondir.folder '\' Sessiondir.name '\*.mat']);
    
    for S=1:length(Sessions)
        load([Sessiondir.folder '\' Sessiondir.name '\' Sessions(S).name])

        %normalize data ********************************************************
        t=linspace(0,length(ECG.values)/500,length(ECG.values));
        tc=linspace(0,length(cECG.values)/500,length(cECG.values));
        ECGdata=(ECG.values - min(ECG.values)) / ( max(ECG.values) - min(ECG.values) );
        cECGdata=(cECG.values - min(cECG.values)) / ( max(cECG.values) - min(cECG.values) );
        ECGdata(isnan(ECGdata)) = []; % Removing nans        
        cECGdata(isnan(cECGdata)) = []; % Removing nans

        %Calc R peaks ********************************************************
        r=length(cECGdata);
        WindowWidth=500*30;

        padding=0; %Determine if the RR should be same length as ECG. Don`t have to be
        plotting=0; %plotting Ralphs RR detection

        RR_idx_cell=cell(ceil(length(cECGdata)/WindowWidth),1);% creating empty cell
        RR_trace_cell=cell(ceil(length(cECGdata)/WindowWidth),1);% creating empty cell

        RR_idx=[];
        RR_trace=[];
        istFlach=[];
        j=1;
        Ralphsfactor=Ralph(P);
%cECG ****************************
%Ralph
    if Methode=='R'
        for i=1:WindowWidth:(numel(cECGdata)-WindowWidth)    % Calculate The RR peaks in windows, otherwise the big Amplitude difference destroys the RR calc
            Datafut=cECGdata(i:i+WindowWidth);
            if range(Datafut)~=0 %If only one value, then Ralphs detector goes mad              
                [RR_idx_cell{j}, ~, ~, ~, ~, RR_trace_cell{j}, ~] = ecg_find_rpeaks(tc,Ralphsfactor*Datafut, 500, 250,plotting,0); %, , , maxrate,plotting,saving   -1* because Ralph optimized for a step s slope, we also have steep Q slope. inverting fixes that probel 
                RR_idx_cell{j}=RR_idx_cell{j}+(j-1)*WindowWidth;
            else 
                RR_idx_cell{j}=nan;
                RR_trace_cell{j}=nan;  
                istFlach=[istFlach,j];
            end
            j=j+1;
        end
        
        for i =1:length(RR_idx_cell) % Stich the Windows together again
            RR_idx=[RR_idx,RR_idx_cell{i}];
        end
       
        for i=1:length(RR_trace_cell)% Stich the Windows together again
            RR_trace=[RR_trace, RR_trace_cell{i}];
        end   
%Michiel
    elseif Methode=='M'
        cECGdata(isnan(cECGdata)) = []; % Removing nans
        
        [RR_idxM] = streamingpeakdetection(cECGdata', 500, [60 256], plotting, 18.5, 1024);
        RR_trace=diff(RR_idxM.peakPositionArray)./500; % Calculating the time between the R peaks in seconds
        RR_idx=RR_idxM.peakPositionArray;      
    end        

%ECG *****************************
%Ralph
    if Methode=='R'
        RR_idx_cell=cell(ceil(length(cECGdata)/WindowWidth),1); % creating empty cell
        RR_trace_cell=cell(ceil(length(cECGdata)/WindowWidth),1);% creating empty cell
        RR_idx2=[];
        RR_trace2=[];
        istFlach2=[];
        for i=1:WindowWidth:(numel(ECGdata)-WindowWidth)     % Calculate The RR peaks in windows, otherwise the big Amplitude difference destroys the RR calc
            Datafut2=ECGdata(i:i+WindowWidth);
            if range(Datafut2)~=0 %If only one value, then Ralphs detector goes mad               
                [RR_idx_cell{j}, ~, ~, ~, ~, RR_trace_cell{j}, ~] = ecg_find_rpeaks(tc,Ralphsfactor*Datafut2, 500, 250,plotting,0); %, , , maxrate,plotting,saving   -1* because Ralph optimized for a step s slope, we also have steep Q slope. inverting fixes that probel 
                RR_idx_cell{j}=RR_idx_cell{j}+(j-1)*WindowWidth;                      
            else
                RR_idx_cell{j}=nan;
                RR_trace_cell{j}=nan;
                istFlach2=[istFlach2,j];
            end
            j=j+1;
        end
% 
%         [~,ia,ib]=intersect(istFlach,istFlach2); % ifboth have the same nan,delete both
%         istFlach(ia)=[];istFlach2(ib)=[];
%         
%         for F=1:length(istFlach2)
%             delete
        

        for i =1:length(RR_idx_cell)% Stich the Windows together again
            RR_idx2=[RR_idx,RR_idx_cell{i}];
        end
        for i=1:length(RR_trace_cell)% Stich the Windows together again
            RR_trace2=[RR_trace2, RR_trace_cell{i}];
        end


%Michiel
    elseif Methode=='M'
        ECGdata(isnan(ECGdata)) = []; % Removing nans
        
        [RR_idxM] = streamingpeakdetection(ECGdata', 500, [60 256], plotting, 18.5, 1024);
        RR_trace2=diff(RR_idxM.peakPositionArray)./500; % Calculating the time between the R peaks in seconds
        RR_idx2=RR_idxM.peakPositionArray;      
    end        
        
        RR_trace(isnan(RR_trace)) = 0; % Removing nans
        RR_trace2(isnan(RR_trace2)) = 0;
        
        %Calc correlation between RR signals of ECG and cECG signals ********************************************************

        [C1,lag1]=xcorr(RR_trace2,RR_trace);
        % C1=C1/(norm(RR_trace2)*norm(RR_trace)); %normalizing to 0-1 for unequal lengthed signals. https://nl.mathworks.com/matlabcentral/answers/5275-algorithm-for-coeff-scaling-of-xcorr
        C1=C1/(sqrt(sum(abs(RR_trace2).^2)*sum(abs(RR_trace).^2)));%normalizing to 0-1 for unequal lengthed signals. 
        [~,idxmax]=find((C1)==max(C1)); % find the inddex where the correlation is highest
        lag=lag1(idxmax);% The actual shift of both signals
        
        Corr_Session_Cell{P,1}(S)=max(C1);
        Corr_Session_Cell{P,2}(S)=lag;
        Corr_Session(S)=max(C1);

    end % for S in Sessions
    Corr_Patient(P,1)=mean(Corr_Session);
    Corr_Patient(P,2)=std(Corr_Session);
    

end % for P in Patients



clearvars -except Corr_Patient Corr_Session_Cell Corr_Session











