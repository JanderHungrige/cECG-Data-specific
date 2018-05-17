function [TPArray, FPArray, FNArray, results] = evaluate_detection_performance_per_motion_level(fid, RRpositions, RRpositionsREF, motionUp, FS, initOffset)
% And write it in the text file fid


[TP, FP, FN, TPArray, FPArray, FNArray, SD, offsetArray] = DeterminePeakOffsets(RRpositions,RRpositionsREF, FS, 0.15, 0);%-round(25e-3*250)

combo=[[zeros(1,length(TPArray));TPArray],[3*ones(1,length(FPArray));FPArray],[ones(1,length(FNArray));FNArray]];%TP=0, FP=3, FN=1
[xxx,index]=sort(combo(2,:),2);
sortedCombo=combo(:,index);
diffCombo=diff(sortedCombo,[],2);
diffCombo(1,:)=diffCombo(1,:)+sortedCombo(1,2:end); %the first line contains zeros where the intervals (2nd line) are reliable

availableHRTime=sum(diffCombo(2,diffCombo(1,:)==0));%in samples
totalTime=sum(diffCombo(2,:)); %RRpositions(end)-RRpositionsREF(1); %the comparison stops at the end of RRpositions
%disp(['For a recording of ' num2str(totalTime./250/60) ' minutes, the iHR was available ' num2str(availableHRTime/totalTime*100) '% of the time.'])
%disp([num2str(round(TP/(FN+TP)*100)) '% of the R-peaks were correctly detected (TP), ' num2str(round(FP/(FN+TP)*100)) '% were added (FP) and ' num2str(round(FN/(FN+TP)*100)) '% were missed (FN).'])
SEN=TP/(TP+FN);
PPV=TP/(TP+FP);
ERROR=(FP + FN) / (TP + FN);
% Length of unreliable periods
[meanLength, maxLength, countLong, saved_ilength]=count_lost_segments(diffCombo, FS);
fprintf(fid,'%s\t\t%3.2f\t\t%3i\t%3i\t\t%3.1f\t%1.2f\t\t\t%2.1f\t%2.1f\t%1.2f\n','Whole data',totalTime/FS/60,round(availableHRTime/totalTime*100), countLong, maxLength, meanLength, PPV*100, SEN*100, ERROR);
fprintf(fid, '\n');
results.wholeData=[totalTime/500/60,availableHRTime/FS/60, countLong, maxLength, meanLength, PPV*100, SEN*100, ERROR, saved_ilength];

% Segment data based on motion
%----------------------------

% Upsample sortedCombo
% sortedComboUp=zeros([2 length(motionUp)]);
% sortedComboUp(1,:)=sortedComboUp(1,:)-1; %first line cannot be full of zeros since 0 means TP
% sortedComboUp(:,sortedCombo(2,:)-STARTtime*60*250)=sortedCombo;

% Downsample motion (keep one motion indicator for each peak )
%CHANGED BY JAN exchange of commented function
motion=motionUp(round(sortedCombo(2,:)-initOffset*FS)); %USEE THIS OR THE ONE BELOW
% motion=motionUp((sortedComboUp(2,:)~=0)); % USE THIS OR THE ONE A LINE UP
%CHanged by jan 
% Segmentation
%----> No Motion:
% peak selection:
% TP=sum(sortedComboUp(1,:)==0 & motionUp==1);
% FP=sum(sortedComboUp(1,:)==3 & motionUp==1);
% FN=sum(sortedComboUp(1,:)==1 & motionUp==1);
TP=sum(sortedCombo(1,:)==0 & motion==1);
FP=sum(sortedCombo(1,:)==3 & motion==1);
FN=sum(sortedCombo(1,:)==1 & motion==1);

if TP>1 || FP>1 || FN>1
    SEN=TP/(TP+FN);
    PPV=TP/(TP+FP);
    ERROR=(FP + FN) / (TP + FN);
    % RR interval selection:
    selectedCombo=diffCombo(:,motion(2:end)==1); % motion(2:end) = associated motion level (taken from the second peak of the interval)
    availableHRTime=sum(selectedCombo(2,selectedCombo(1,:)==0));%in samples
    totalTime=sum(selectedCombo(2,:)); %the comparison stops at the end of RRpositions
    % Length of unreliable periods
    [meanLength, maxLength, countLong, saved_ilength]=count_lost_segments(selectedCombo, FS);
    fprintf(fid,'%s\t\t%3.2f\t\t%3i\t%3i\t\t%3.1f\t%1.2f\t\t\t%2.1f\t%2.1f\t%1.2f\n','No motion',totalTime/FS/60,round(availableHRTime/totalTime*100), countLong, maxLength, meanLength, PPV*100, SEN*100, ERROR);
    results.noMotion=[totalTime/250/60,availableHRTime/FS/60, countLong, maxLength, meanLength, PPV*100, SEN*100, ERROR, saved_ilength];
else
    results.noMotion=zeros(1,8);
end

%----> Low Motion:
% peak selection:
TP=sum(sortedCombo(1,:)==0 & motion==2);
FP=sum(sortedCombo(1,:)==3 & motion==2);
FN=sum(sortedCombo(1,:)==1 & motion==2);
if TP>1 || FP>1 || FN>1
    SEN=TP/(TP+FN);
    PPV=TP/(TP+FP);
    ERROR=(FP + FN) / (TP + FN);
    % RR interval selection:
    selectedCombo=diffCombo(:,motion(2:end)==2); % motion(2:end) = associated motion level (taken from the second peak of the interval)
    availableHRTime=sum(selectedCombo(2,selectedCombo(1,:)==0));%in samples
    totalTime=sum(selectedCombo(2,:)); %the comparison stops at the end of RRpositions
    % Length of unreliable periods
    [meanLength, maxLength, countLong, saved_ilength]=count_lost_segments(selectedCombo, FS);
    fprintf(fid,'%s\t\t%3.2f\t\t%3i\t%3i\t\t%3.1f\t%1.2f\t\t\t%2.1f\t%2.1f\t%1.2f\n','Low motion',totalTime/FS/60,round(availableHRTime/totalTime*100), countLong, maxLength, meanLength, PPV*100, SEN*100, ERROR);
    results.lowMotion=[totalTime/250/60,availableHRTime/FS/60, countLong, maxLength, meanLength, PPV*100, SEN*100, ERROR, saved_ilength];
else
    results.lowMotion=zeros(1,8);

end

%----> High Motion:
% peak selection:
TP=sum(sortedCombo(1,:)==0 & motion==3);
FP=sum(sortedCombo(1,:)==3 & motion==3);
FN=sum(sortedCombo(1,:)==1 & motion==3);
if TP>1 || FP>1 || FN>1
    SEN=TP/(TP+FN);
    PPV=TP/(TP+FP);
    ERROR=(FP + FN) / (TP + FN);
    % RR interval selection:
    selectedCombo=diffCombo(:,motion(2:end)==3); % motion(2:end) = associated motion level (taken from the second peak of the interval)
    availableHRTime=sum(selectedCombo(2,selectedCombo(1,:)==0));%in samples
    totalTime=sum(selectedCombo(2,:)); %the comparison stops at the end of RRpositions
    % Length of unreliable periods
    [meanLength, maxLength, countLong, saved_ilength]=count_lost_segments(selectedCombo, FS);
    fprintf(fid,'%s\t\t%3.2f\t\t%3i\t%3i\t\t%3.1f\t%1.2f\t\t\t%2.1f\t%2.1f\t%1.2f\n','High motion',totalTime/FS/60,round(availableHRTime/totalTime*100), countLong, maxLength, meanLength, PPV*100, SEN*100, ERROR);
    fprintf(fid, '\n');
    results.highMotion=[totalTime/250/60,availableHRTime/FS/60, countLong, maxLength, meanLength, PPV*100, SEN*100, ERROR, saved_ilength];
else
    results.highMotion=zeros(1,8);
end

%----> Intervention:
% peak selection:
TP=sum(sortedCombo(1,:)==0 & (motion==4 | motion==5));
FP=sum(sortedCombo(1,:)==3 & (motion==4 | motion==5));
FN=sum(sortedCombo(1,:)==1 & (motion==4 | motion==5));
if TP>1 || FP>1 || FN>1
    SEN=TP/(TP+FN);
    PPV=TP/(TP+FP);
    ERROR=(FP + FN) / (TP + FN);
    % RR interval selection:
    selectedCombo=diffCombo(:,(motion(2:end)==4|motion(2:end)==5)); % motion(2:end) = associated motion level (taken from the second peak of the interval)
    availableHRTime=sum(selectedCombo(2,selectedCombo(1,:)==0));%in samples
    totalTime=sum(selectedCombo(2,:)); %the comparison stops at the end of RRpositions
    % Length of unreliable periods
    [meanLength, maxLength, countLong, saved_ilength]=count_lost_segments(selectedCombo, FS);
    fprintf(fid,'%s\t\t%3.2f\t\t%3i\t%3i\t\t%3.1f\t%1.2f\t\t\t%2.1f\t%2.1f\t%1.2f\n','Intervention',totalTime/FS/60,round(availableHRTime/totalTime*100), countLong, maxLength, meanLength, PPV*100, SEN*100, ERROR);
    fprintf(fid, '\n');
    results.intervention=[totalTime/250/60,availableHRTime/FS/60, countLong, maxLength, meanLength, PPV*100, SEN*100, ERROR, saved_ilength];
else
    results.intervention=zeros(1,8);
    
end

%----> Ghost:
% peak selection:
TP=sum(sortedCombo(1,:)==0 & motion==6);
FP=sum(sortedCombo(1,:)==3 & motion==6);
FN=sum(sortedCombo(1,:)==1 & motion==6);
if TP>1 || FP>1 || FN>1
    SEN=TP/(TP+FN);
    PPV=TP/(TP+FP);
    ERROR=(FP + FN) / (TP + FN);
    % RR interval selection:
    selectedCombo=diffCombo(:,motion(2:end)==6); % motion(2:end) = associated motion level (taken from the second peak of the interval)
    availableHRTime=sum(selectedCombo(2,selectedCombo(1,:)==0));%in samples
    totalTime=sum(selectedCombo(2,:)); %the comparison stops at the end of RRpositions
    % Length of unreliable periods
    [meanLength, maxLength, countLong, saved_ilength]=count_lost_segments(selectedCombo, FS);
    fprintf(fid,'%s\t\t%3.2f\t\t%3i\t%3i\t\t%3.1f\t%1.2f\t\t\t%2.1f\t%2.1f\t%1.2f\n','Ghosts',totalTime/FS/60,round(availableHRTime/totalTime*100),countLong, maxLength, meanLength, PPV*100, SEN*100, ERROR);
    fprintf(fid, '\n');
    results.ghosts=[totalTime/250/60,availableHRTime/FS/60, countLong, maxLength, meanLength, PPV*100, SEN*100, ERROR, saved_ilength];
else
    results.ghosts=zeros(1,8);
end

%----> Ignored data:
% peak selection:
TP=sum(sortedCombo(1,:)==0 & motion==7);
FP=sum(sortedCombo(1,:)==3 & motion==7);
FN=sum(sortedCombo(1,:)==1 & motion==7);
if TP>1 || FP>1 || FN>1
    SEN=TP/(TP+FN);
    PPV=TP/(TP+FP);
    ERROR=(FP + FN) / (TP + FN);
    % RR interval selection:
    selectedCombo=diffCombo(:,motion(2:end)==7); % motion(2:end) = associated motion level (taken from the second peak of the interval)
    availableHRTime=sum(selectedCombo(2,selectedCombo(1,:)==0));%in samples
    totalTime=sum(selectedCombo(2,:)); %the comparison stops at the end of RRpositions
    % Length of unreliable periods
    [meanLength, maxLength, countLong, saved_ilength]=count_lost_segments(selectedCombo, FS);
    fprintf(fid,'%s\t\t%3.2f\t\t%3i\t%3i\t\t%3.1f\t%1.2f\t\t\t%2.1f\t%2.1f\t%1.2f\n','Ignored',totalTime/FS/60,round(availableHRTime/totalTime*100), countLong, maxLength, meanLength, PPV*100, SEN*100, ERROR);
    fprintf(fid, '\n');
    results.ignored=[totalTime/250/60,availableHRTime/FS/60, countLong, maxLength, meanLength, PPV*100, SEN*100, ERROR, saved_ilength];
else
    results.ignored=zeros(1,8);
end

