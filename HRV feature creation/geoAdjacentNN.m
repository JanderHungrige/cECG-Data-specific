% geometrical analysis 2
% histogram
% sample density distribution of adjacent NN interval durations
%HRVindex
%TINN

function geoAdjacentNN(pat,saving,win,plot)
% pat=[3,4,5,6,7,9,10,13,15];
pat=13;
for i=1:length(pat)
    Neonate=pat(i);
    
%AS
    load(['C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\RR distance2\RtimeAS' num2str(Neonate)]);
    if timebetweenR
        %remove outlier
        all_idx = 1:length(timebetweenR);
        outlier_idx = abs(timebetweenR - median(timebetweenR)) > 10*std(timebetweenR); %| abs(y - median(y)) > 3*std(y) % Find outlier idx
        timebetweenR(outlier_idx) = interp1(all_idx(~outlier_idx), timebetweenR(~outlier_idx), all_idx(outlier_idx)); % Linearly interpolate over outlier idx for x
        
        %calcuate difference of adjacent NN intervals
        NNadjcacent= diff(timebetweenR);
        
         bin=floor((max(NNadjcacent)-min(NNadjcacent))/(8/1000));   %bins for geometrical histogram analysis should be optimal 8ms (7,8125ms) broad (adult value)
         [histo,bin]=hist(NNadjcacent,bin);  %creating histogram with n bins
         histo=(histo./sum(histo)).*100; % normalizing the histogram in %
    end
%QS
    load(['C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\RR distance2\RtimeQS' num2str(Neonate)]);
    if timebetweenR

        %remove outlier
        all_idx = 1:length(timebetweenR);
        outlier_idx = abs(timebetweenR - median(timebetweenR)) > 10*std(timebetweenR); %| abs(y - median(y)) > 3*std(y) % Find outlier idx
        timebetweenR(outlier_idx) = interp1(all_idx(~outlier_idx), timebetweenR(~outlier_idx), all_idx(outlier_idx)); % Linearly interpolate over outlier idx for x
  
        %calcuate difference of adjacent NN intervals
        NNadjcacent= diff(timebetweenR);
        
       %  bin2=floor((max(NNadjcacent)-min(NNadjcacent))/(8/1000));   %bins for geometrical histogram analysis should be optimal 8ms (7,8125ms) broad (adult value)
        [histo2]=hist(NNadjcacent,bin);
        histo2=(histo2./sum(histo2)).*100;

        %statistical comparison between AS&QS histogram with p value
        [h,p(i)] = kstest2(histo,histo2);
    end
    
figure
hold on
bar(bin,[histo', histo2'])
% bar(bin,histo,'b')
% bar(bin2,histo2,'y')
text(10,37,4,['pvalue:' num2str(p)])
xlabel('RR distance in S');
ylabel('Appearance in percent');
legend('AS', 'QS')
title(['Patient' num2str(Neonate)])
xlim([-0.1 0.1])
hold off
end