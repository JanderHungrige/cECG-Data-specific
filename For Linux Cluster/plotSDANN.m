% %SDANN
function plotSDANN(plotting, win, pat, overlap)
if plotting
    for j=1:length(win)
        for o=1:length(overlap)
            a=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];%
% a=[3 4 5 6 7 9]
                cd('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\windowed_ECG_chunks\HRV analysis\Ralphs')

%%%%%%%%%%%%%%%%AS
                if exist(fullfile(cd, ['SDANN_AS_' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,o)) '.mat']), 'file') == 2 % ==> 0 or 2
                   load((fullfile(cd, ['SDANN_AS_' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,o)) '.mat'])))

                   SDANNAS(SDANNAS==0)=nan; % changing zeroes to nan. Nans are not plotted
%                    SDANNAScollect(j,:)=SDANNAS;

                    figure(2)
                    plot(a,SDANNAS,'x b')
%                     title(['Patient ' num2str(pat(1,i)) ' with Window size: ' num2str(win(1,j)) ' overlap ' num2str(overlap(1,o))])
                    hold on %continued in QS
                
%bar plot of histogram
%                    bin=20;
%                    figure(2)
%                         bar(edges(1,1:end-1),histoc)
%                         hold on 
%                        
%                         [histoc,edges]=histcounts(SDANNAScollect,bin);
%                         histoc=(histoc./sum(histoc)).*100; % normalizing the histogram in %
                end

%%%%%%%%%%%%%%%QS
                if exist(fullfile(cd, ['SDANN_QS_' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,o)) '.mat']), 'file') == 2 % ==> 0 or 2
                   load((fullfile(cd, ['SDANN_QS_' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,o)) '.mat'])))

                   SDANNQS(SDANNQS==0)=nan;
%                    SDANNQScollect(j,:)=SDANNQS;

                   SDANNQS(1,14:15)=nan; % QS is missing in 15. TO plot it it has to be added manually
                    plot(a,SDANNQS,'x r')
                    
                    ylabel('SDANN');
                    xlabel('Patient')
                   % title(['win ' num2str(win(1,j)) ' overlap ' num2str(overlap(1,o))]);
                    legend('AS', 'QS')
                    xlim([2 16])
                    x=[ 3 4 5 6 7 9 13 15]; % new x axis
                    set(gca,'XTick',x); % Change x-axis ticks
                       
% bar plot of histogram                        
%                         [histoc2,edges2]=histcounts(SDANNQScollect,bin);
%                          histoc2=(histoc2./sum(histoc2)).*100; % normalizing the histogram in %                                                               
%                         
%                         bar(edges2(1,1:end-1),histoc2,'y')
%                         hold off
%                         
%                         title(['Distribution of SDANNS'])
%                         xlabel('SDANN')
%                         ylabel('Percent')
%                         xlim([0 0.25])
                end
            end
        end
    end
end
