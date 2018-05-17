function plotRMSSD(plotting, win, pat, overlap)
if plotting
    for j=1:length(win)
        for o=1%:length(overlap)
                cd('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\windowed_ECG_chunks\HRV analysis\Ralphs\')
                if exist(fullfile(cd, ['RMSSD' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,o)) '.mat']), 'file') == 2 % ==> 0 or 2
                   load((fullfile(cd, ['RMSSD' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,o)) '.mat'])))

RMSSD_AS(cellfun(@isempty, RMSSD_AS)) = {NaN};
RMSSD_QS(cellfun(@isempty, RMSSD_QS)) = {NaN};
x1=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
x2=[1.3 2.3 3.3 4.3 5.3 6.3 7.3 8.3 9.3 10.3 11.3 12.3 13.3 14.3 15.3];
                   RMSSD_AS=cell2mat(RMSSD_AS);
                   RMSSD_QS=cell2mat(RMSSD_QS);
                   RMSSD_QS(14:15,:)=nan; %QS is shorter that AS therfore, manually nans have to be filled in
figure(2)
                   plot(x1,RMSSD_AS./1000,'x b')
                   hold on 
                   plot(x2,RMSSD_QS./1000, 'x r')
                   ylim([0 .500])
                   x=[3 4 5 6 7 9 13 15];
                   set(gca,'XTick',x); % Change x-axis ticks

                       title(['Distribution of RMSSD'])
                        ylabel('RMSSD [s]')
                        xlabel('Patient')
                        
                end
        end
    end
    
end