% geometrical analysis
% histogram
% sample density distribution of NN interval durations 
%HRVindex
%TINN


function geoNN(pat,saving,win,plot)
    for i=1:length(pat)
        Neonate=pat(i);

    %AS
        load(['C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\RR distance2\RtimeAS' num2str(Neonate)]);
        if timebetweenR
            %remove outlier
            all_idx = 1:length(timebetweenR);
            outlier_idx = abs(timebetweenR - median(timebetweenR)) > 10*std(timebetweenR); %| abs(y - median(y)) > 3*std(y) % Find outlier idx
            timebetweenR(outlier_idx) = interp1(all_idx(~outlier_idx), timebetweenR(~outlier_idx), all_idx(outlier_idx)); % Linearly interpolate over outlier idx for x

             bin=floor((max(timebetweenR)-min(timebetweenR))/(8/1000));   %bins for geometrical histogram analysis should be optimal 8ms (7,8125ms) broad (adult value)
             [histo,bin]=hist(timebetweenR,bin);  %creating histogram with 100 bins

            %calculting HRV triangular index
            Y_AS=max(histo);
            Triangularidx_AS(i,1)=(length(timebetweenR)/Y_AS);

            %calculating TINN
            Y=normpdf(histo);
            YY=mvnpdf(histo);
       histfit(timebetweenR,134)
             figure
            bar(bin,histo)
          %   histo=(histo./sum(histo)).*100; % normalizing the histogram in %
        end
    %     pd=fitdist(histo','Normal');
    %     yval=pdf(pd,bin);

    %QS
        load(['C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\RR distance2\RtimeQS' num2str(Neonate)]);
        if timebetweenR

            %remove outlier
            all_idx = 1:length(timebetweenR);
            outlier_idx = abs(timebetweenR - median(timebetweenR)) > 10*std(timebetweenR); %| abs(y - median(y)) > 3*std(y) % Find outlier idx
            timebetweenR(outlier_idx) = interp1(all_idx(~outlier_idx), timebetweenR(~outlier_idx), all_idx(outlier_idx)); % Linearly interpolate over outlier idx for x

           % bin2=floor((max(timebetweenR)-min(timebetweenR))/(8/1000));   %bins for geometrical histogram analysis should be optimal 8ms (7,8125ms) broad (adult value)
            [histo2]=hist(timebetweenR,bin);
            histo2=(histo2./sum(histo2)).*100;

            %statistical comparison between AS&QS histogram with p value
            [h,p(i)] = kstest2(histo,histo2);
        end

        %calculting HRV triangular index
        Y_QS=max(histo2);
        Triangularidx_QS(i,1)=(length(timebetweenR)/Y_QS);


        plot=0; % for figure plot=1
     if plot  
    figure
    hold on
    bar(bin,[histo',histo2'])
    % bar(bin,histo,'b')
    % bar(bin2,histo2,'y')
    text(10,37,4,['pvalue:' num2str(p)])
    xlabel('RR distance in S');
    ylabel('Appearance in percent');
    legend('AS', 'QS')
    title(['Patient' num2str(Neonate)])
    % xlim([0.3 0.5])
    hold off
     end
    end
end