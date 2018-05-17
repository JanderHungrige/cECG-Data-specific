function [MeanArclength, SDarclength]=SDarclength(Neonate,pat,arclenAS,arclenQS)
if Neonate==15
    temp=(reshape(arclenAS',1,[]));
    temp2=(reshape(arclenQS',1,[]));

    for j=1:(length(arclenAS))-10
    meanAclenAS(j,1)=nanmean(temp(1,j:j+9));
    meanAclenAS(j,2)= nanstd(temp(1,j:j+9));
    end
    
    for j=1:(length(arclenQS))-10
    meanAclenQS(j,1)=nanmean(temp2(1,j:j+9));
    meanAclenQS(j,2)= nanstd(temp2(1,j:j+9));
    end

 
    
%%%% HISTOGRAM    
        
[counts, binValue]= hist(meanAclenAS(:,1),200); %D is your data and 140 is number of bins.
counts = 100*counts/sum(counts); % normalize to unit length. Sum of h now will be 1.
[counts2, binValue2]= hist(meanAclenQS(:,1),200); %D is your data and 140 is number of bins.
counts2 = 100*counts2/sum(counts2); % normalize to unit length. Sum of h now will be 1.

   bar(binValue, counts); hold on 
%    xlim([0 1000])
   bar(binValue2, counts2, 'y'); hold off
%    xlim([0 1000])
   title('Mean arclength')

[counts, binValue]= hist(meanAclenAS(:,2),200); %D is your data and 140 is number of bins.
counts = 100*counts/sum(counts); % normalize to unit length. Sum of h now will be 1.
[counts2, binValue2]= hist(meanAclenQS(:,2),200); %D is your data and 140 is number of bins.
counts2 = 100*counts2/sum(counts2); % normalize to unit length. Sum of h now will be 1.

   bar(binValue, counts); hold on 
%    xlim([0 1000])
   bar(binValue2, counts2, 'y'); hold off
%    xlim([0 1000])
title('std linelength')
end %if Neonate==15
end