
%Drag and drop the DAQ files from E:\cECG_study\For_Quick_Annotator\participant x ... into matlab first
plotr=1 % plotting the r peaks


%normalize data ********************************************************
t=linspace(0,length(ECG.values)/500,length(ECG.values));
tc=linspace(0,length(cECG.values)/500,length(cECG.values));
ECGdata=(ECG.values - min(ECG.values)) / ( max(ECG.values) - min(ECG.values) );
cECGdata=(cECG.values - min(cECG.values)) / ( max(cECG.values) - min(cECG.values) );

%Calc R peaks ********************************************************
r=length(cECGdata);
WindowWidth=500*10;
Ralphsfactor={1;-1;-1;1;1;-1;1;-1; 1; 1;-1; 1;-1; 1;-1;-1;-1;-1};%Determine if the ECG signal should be turned -1 or not 1. 
             %1  2 3  4 5  6 7  8  9  10 11 12 13 14 15 16 17 18
Ralphsfactor=1;
j=1;
if plotr
    padding=0; %Determine if the RR should be same length as ECG. Don`t have to be
    plotting=0; %plotting Ralphs RR detection

    RR_idx_cell=cell(ceil(length(cECGdata)/WindowWidth),1);% creating empty cell
    RR_trace_cell=cell(ceil(length(cECGdata)/WindowWidth),1);% creating empty cell

    RR_idx=[];
    RR_trace=[];
    for i=1:WindowWidth:(numel(cECGdata)-WindowWidth+1)  % Calculate The RR peaks in windows, otherwise the big Amplitude difference destroys the RR calc  
        Datafut=cECGdata(i:i+WindowWidth);
        [RR_idx_cell{j}, ~, ~, ~, ~, RR_trace_cell{j}, ~] = ecg_find_rpeaks(tc,Ralphsfactor*Datafut, 500, 250,plotting,0); %, , , maxrate,plotting,saving   -1* because Ralph optimized for a step s slope, we also have steep Q slope. inverting fixes that probel 
        RR_idx_cell{j}=RR_idx_cell{j}+(j-1)*WindowWidth;
        j=j+1;
    end
    for i =1:length(RR_idx_cell)% Stich the Windows together again
        RR_idx=[RR_idx,RR_idx_cell{i}];
    end
    for i=1:length(RR_trace_cell)% Stich the Windows together again
        RR_trace=[RR_trace, RR_trace_cell{i}];
    end
    R=NaN(length(cECGdata),1);
    R(RR_idx)=cECGdata(RR_idx)+0.6; % creating R signal to lay over ECG signal for plotting
    
    
    
    RR_idx_cell=cell(ceil(length(cECGdata)/WindowWidth),1);% creating empty cell
    RR_trace_cell=cell(ceil(length(cECGdata)/WindowWidth),1);% creating empty cell

    RR_idx2=[];
    RR_trace2=[];
    for i=1:WindowWidth:(numel(ECGdata)-WindowWidth+1)    % Calculate The RR peaks in windows, otherwise the big Amplitude difference destroys the RR calc
        Datafut=ECGdata(i:i+WindowWidth);
        [RR_idx_cell{j}, ~, ~, ~, ~, RR_trace_cell{j}, ~] = ecg_find_rpeaks(tc,Ralphsfactor*Datafut, 500, 250,plotting,0); %, , , maxrate,plotting,saving   -1* because Ralph optimized for a step s slope, we also have steep Q slope. inverting fixes that probel 
        RR_idx_cell{j}=RR_idx_cell{j}+(j-1)*WindowWidth;
        j=j+1;
    end
    for i =1:length(RR_idx_cell)% Stich the Windows together again
        RR_idx2=[RR_idx,RR_idx_cell{i}];
    end
    for i=1:length(RR_trace_cell)% Stich the Windows together again
        RR_trace2=[RR_trace2, RR_trace_cell{i}];
    end
    R2=NaN(length(ECGdata),1);
    R2(RR_idx)=ECGdata(RR_idx); % creating R signal to lay over ECG signal for plotting

end
    

RR_trace(isnan(RR_trace)) = 0;% Removing nans
RR_trace2(isnan(RR_trace2)) = 0;

%Calc correlation between RR signals of ECG and cECG signals ********************************************************

[C1,lag1]=xcorr(RR_trace2,RR_trace);
% C1=C1/(norm(RR_trace2)*norm(RR_trace)); %normalizing to 0-1 for unequal
% lengthed signals. https://nl.mathworks.com/matlabcentral/answers/5275-algorithm-for-coeff-scaling-of-xcorr
C1=C1/(sqrt(sum(abs(RR_trace2).^2)*sum(abs(RR_trace).^2)));
% C2=corrcoef(RR_trace,RR_trace2);
figure;plot(lag1/500,C1,'g')

% Plot ********************************************************
figure
set(gcf,'color','w')

% a=subplot(2,1,1);
hECG=line(t,ECGdata);
set(gca,'yticklabel','')
% b=subplot(2,1,2);
hcECG=line(tc,cECGdata+0.5);
if plotr
    hR=line(tc,R);
end
% linkaxes([a,b],'x')

% set(a,axes,'Color','none','XColor','none' )
if plotr
    set(hR                            , ...
      'LineStyle'       , ': '      , ...
      'Color'           , [1 0 0]  ,...
      'Marker'       , 'x');
end

set(hcECG                            , ...
  'LineStyle'       , ': '      , ...
  'Color'           , [0.4 0.4 0.4]  ,...
  'LineWidth'       , 1.5         );

set(hECG                            , ...
  'LineStyle'       , '-'      , ...
  'Color'           , [0 0 0]  ,...
  'LineWidth'       , 1           );

htitle=title('');
hxlabel=xlabel('time[s]'                     );
% hylabel=ylabel('Amplitude normalized'                      );
hLegend = legend([hECG,hcECG,hR],...
    'ECG',...
    'cECG',...
    'Deteted cECG R peaks');

set(hLegend, ...
    'FontSize'   , 8          , ...
    'FontName'   , 'Helvetica');

set([htitle, hxlabel],...%, hylabel], ...
    'FontSize'   , 12          , ...
    'FontName'   , 'AvantGarde');

set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 finalPlot1.eps




















