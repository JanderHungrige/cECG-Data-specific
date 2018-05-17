% %SDNN
function plotSDNN(plotting, win, pat, overlap)
if plotting
    for j=1:length(win)
        for o=1%:length(overlap)
                cd('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\windowed_ECG_chunks\HRV analysis')
                if exist(fullfile(cd, ['SDNN_' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,o)) '.mat']), 'file') == 2 % ==> 0 or 2
            
                    load((fullfile(cd, ['SDNN_' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,o)) '.mat'])))
                end
%plotting dsitribution                 
                    for i=1:length(pat)
%                         %  histo=hist(cell2mat(SDNNAS(pat(i),:)),100);
%                        % histo=(histo./sum(histo)).*100; % normalizing the histogram in %
%   
%                        bin=150;
% 
%                                              
%                         [histoc,edges]=histcounts(cell2mat(SDNNAS(pat(i),:)),bin);
%                         histoc=(histoc./sum(histoc)).*100; % normalizing the histogram in %
%                         
%                         [histoc2,edges2]=histcounts(cell2mat(SDNNQS(pat(i),:)),bin);
%                         histoc2=(histoc2./sum(histoc2)).*100; % normalizing the histogram in %
%                         
%                         figure
% %                         bar([edges edges2],[histoc', histoc2'])
%                         bar(edges(1,1:end-1),histoc)
%                         hold on 
%                         bar(edges2(1,1:end-1),histoc2,'y')
%                         hold off
%                         
%                         title(['Patient:' num2str(pat(i)) ' window:' num2str(win(1,j)) ' overlap:' num2str(overlap(1,o))])
%                         xlabel('std')
%                         ylabel('percent')
%                         
% %                              bar(histo)
% %                                 histo=(histo./sum(histo)).*100; % normalizing the histogram in %
% %                                 bar(histo)
% figure(1)
%                     plot(cell2mat(SDNNAS(pat(i),:)),'x b')
%                    % title(['Patient ' num2str(pat(1,i)) ' with Window size: ' num2str(win(1,j)) ' overlap ' num2str(overlap(1,o))])
%                    title('SDNN over all patients') 
%                    hold on
%                     plot(cell2mat(SDNNQS(pat(i),:)),'x r')
%                    % hold off
%                     ylabel('SDNN');
%                     xlabel('index')
%                     legend('AS', 'QS')
%                     ylim([0 0.05])
% % 
% %  
% % 
% %                  figure(i)
% %                     subplot(3,2,j)
% %                     plot(cell2mat(SDNNAS(pat(i),:)),'x b')
% %                     title(['Patient ' num2str(pat(1,i)) ' with Window size: ' num2str(win(1,j)) ' overlap ' num2str(overlap(1,o))])
% %                    hold on
% %                     plot(cell2mat(SDNNQS(pat(i),:)),'x r')
% %                     hold off
% %                     ylabel('SDNN');
% %                     xlabel('index')
% %                     legend('AS', 'QS')
% %                     ylim([0 0.05])
% 
% 
% % 
% %          bin=100;   %bins for geometrical histogram analysis should be optimal 8ms (7,8125ms) broad (adult value)
% %          [histo,bin]=hist(cell2mat(SDNNAS(pat(i),:),bin));  %creating histogram with n bins
% %          histo=(histo./sum(histo)).*100; % normalizing the histogram in %
% %     
% % %QS
% %         [histo2]=hist(cell2mat(SDNNQS(pat(i),:),bin));
% %          histo2=(histo2./sum(histo2)).*100;
% 
%         %statistical comparison between AS&QS histogram with p value
% %         [h,p(i)] = kstest2(histo,histo2);
%     
% %     
% % figure
% % 
% % bar(bin,[histo', histo2'])

%  figure(2)
%         if j==2
%             n=1;
%         subplot(221)
%           plot(cell2mat(SDNNAS(pat(1),:)),'x b')
%                    title(['SDNN pat' num2str(pat(1)) '; win ' num2str(win(1,j))]) 
%                    hold on
%                     plot(cell2mat(SDNNQS(pat(1),:)),'x r')
%                    % hold off
%                     ylabel('SDNN');
%                     xlabel('index')
%                     legend('AS', 'QS')
%                     ylim([0 0.05])
%         end
%                 if j==3
%         subplot(222)
%           plot(cell2mat(SDNNAS(pat(1),:)),'x b')
%                    title(['SDNN all pat' num2str(pat(1)) '; win ' num2str(win(1,j))]) 
%                    hold on
%                     plot(cell2mat(SDNNQS(pat(1),:)),'x r')
%                    % hold off
%                     ylabel('SDNN');
%                     xlabel('index')
%                     legend('AS', 'QS')
%                     ylim([0 0.05])
%                 end
%                   if j==4
%         subplot(223)
%           plot(cell2mat(SDNNAS(pat(1),:)),'x b')
%                    title(['SDNN all pat' num2str(pat(1)) '; win ' num2str(win(1,j))]) 
%                    hold on
%                     plot(cell2mat(SDNNQS(pat(1),:)),'x r')
%                    % hold off
%                     ylabel('SDNN');
%                     xlabel('index')
%                     legend('AS', 'QS')
%                     ylim([0 0.05])
%                   end
%                   if j==5
%         subplot(224)
%           plot(cell2mat(SDNNAS(pat(1),:)),'x b')
%                    title(['SDNN all pat' num2str(pat(1)) '; win ' num2str(win(1,j))]) 
%                    hold on
%                     plot(cell2mat(SDNNQS(pat(1),:)),'x r')
%                    % hold off
%                     ylabel('SDNN');
%                     xlabel('index')
%                     legend('AS', 'QS')
%                     ylim([0 0.05])
%                   end

        end
 figure(3)
                   plot(pat(1),cell2mat(SDNNAS(pat(1),:))./1000,'x b');
                   hold on 
                   plot(pat(1)+0.3,cell2mat(SDNNQS(pat(1),:))./1000,'x r');
                   plot(pat(2),cell2mat(SDNNAS(pat(2),:))./1000,'x b')
                   plot(pat(2)+0.3,cell2mat(SDNNQS(pat(2),:))./1000,'x r')
                   plot(pat(3),cell2mat(SDNNAS(pat(3),:))./1000,'x b')
                   plot(pat(3)+0.3,cell2mat(SDNNQS(pat(3),:))./1000,'x r')
                   plot(pat(4),cell2mat(SDNNAS(pat(4),:))./1000,'x b')
                   plot(pat(4)+0.3,cell2mat(SDNNQS(pat(4),:))./1000,'x r')
                   plot(pat(5),cell2mat(SDNNAS(pat(5),:))./1000,'x b')
                   plot(pat(5)+0.3,cell2mat(SDNNQS(pat(5),:))./1000,'x r')
                   plot(pat(6),cell2mat(SDNNAS(pat(6),:))./1000,'x b')
                   plot(pat(6)+0.3,cell2mat(SDNNQS(pat(6),:))./1000,'x r')
                   ylim([0 1.2])
                   xlim([2 14])
%                    x=[2  4  6  8 10 12 14];
                   x=[ 3 4 5 6 7 8 9 10 11 12 13 ];

                   set(gca,'XTick',x); % Change x-axis ticks
                   title(['SDNN'])
                   ylabel('SDNN')
                   xlabel('Patients')
                   legend('AS', 'QS');

%                figure(2)
%                    plot(NN30_AS,'x b')
%                    hold on 
%                    plot(NN30_QS,'x r')
%                    %ylim([0 0.5])
%                    %x=[3 4 6 7 9 13];
%                    %set(gca,'XTick',x); % Change x-axis ticks
%                    title(['NN30'])
%                    ylabel('NN30')
%                    xlabel('index')
%                figure(3)
%                    plot(NN20_AS,'x b')
%                    hold on 
%                    plot(NN20_QS,'x r')
%                    
%                    %ylim([0 0.5])
%                    %x=[3 4 6 7 9 13];
%                    %set(gca,'XTick',x); % Change x-axis ticks
%                    title(['NN20'])
%                    ylabel('NN20')
%                    xlabel('index')
%                figure(4)
%                    plot(NN10_AS,'x b')
%                    hold on 
%                    plot(NN10_QS,'x r')
%                    
%                    %ylim([0 0.5])
%                    %x=[3 4 6 7 9 13];
%                    %set(gca,'XTick',x); % Change x-axis ticks
%                    title(['NN10'])
%                    ylabel('NN10')
%                    xlabel('index')

                    
              end
         end
    end
end