function [meanLength, maxLength, countLong, saved_ilength]=count_lost_segments(diffCombo, FS)
% diffCombo: 1 st line 0 (correct RRinterval) or 1 (wrong RRinterval)
% 2nd line: RR intevral value (in samples)
% meanLength and maxLength in seconds

tempLength=0;
cumSum=0; numCum=0;
maxLength=0;
countLong=0;
saved_ilength=[];
for i=1:size(diffCombo,2)
    if diffCombo(1,i)==1 %unreliable RR interval
        tempLength=tempLength+diffCombo(2,i); %add the length of the current RR interval
    else
        %update max length
        if tempLength>maxLength
            maxLength=tempLength;
        end
        if tempLength>=10*FS %10 seconds
            countLong=countLong+1;
            saved_ilength=[saved_ilength, tempLength];
            %update mean length
            cumSum=cumSum+tempLength;
            numCum=numCum+1;
        end
        
        %reset length
        tempLength=0;
    end
end
%last one: update max length
if tempLength>maxLength
    maxLength=tempLength;
end
if tempLength>=10*FS %10 seconds
    countLong=countLong+1;
    saved_ilength=[saved_ilength, tempLength];
end
meanLength=cumSum/numCum/FS;
maxLength=maxLength/FS;
saved_ilength=saved_ilength./FS; %in seconds