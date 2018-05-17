
pat=[3,4,5,6,7,9,10,13];
win=[30,60,120,180,240,300];
 %pat=3;
%%SDANN
% % for j=1:length(win)
% %       figure  
% %     for i=1:length(pat)
% %               load(['C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\windowed_ECG_chunks\SDANN_' num2str(win(1,j)) 's_win.mat'])
% %           
% %     plot(cell2mat(SDNNAS(pat(i),:)),'x b')
% %     title(['Patient ' num2str(pat(1,i)) ' with Window size: ' num2str(win(1,j))])
% %     hold on
% %     plot(cell2mat(SDNNQS(pat(i),:)),'x r')
% %     
% %     ylabel('SDANN');
% %     xlabel('index')
% %     legend('AS', 'QS')
% %      ylim([0 0.05])
% %   
% % 
% % 
% % 
% %  figure(i)
% %     subplot(3,2,j)
% %     plot(cell2mat(SDNNAS(pat(i),:)),'x b')
% %     title(['Patient ' num2str(pat(1,i)) ' with Window size: ' num2str(win(1,j))])
% %    hold on
% %     plot(cell2mat(SDNNQS(pat(i),:)),'x r')
% %     hold off
% %     ylabel('SDNN');
% %     xlabel('index')
% %     legend('AS', 'QS')
% %     ylim([0 0.05])
% % 
% % 
% %     end
% %      title([ ' with Window size: ' num2str(win(1,j))])
% % end
% % 

%RMSSD
% for j=1:length(win)
%       figure  
%     for i=1:length(pat)
%               load(['C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\windowed_ECG_chunks\RMSSD' num2str(win(1,j)) 's_win.mat'])
%           
%    fiure(i)
%     plot(cell2mat(RMSSD_QS(pat(i),:)),'x b')
%     title(['Patient ' num2str(pat(1,i)) ' with Window size: ' num2str(win(1,j))])
%     hold on
%     plot(cell2mat(RMSSD_AS(pat(i),:)),'x r')
%     hold off
%     ylabel('RMSSD');
%     xlabel('index')
%     legend('AS', 'QS')
%      ylim([0 0.05])
%   
% 
% 
% 
%  figure(i)
%     subplot(3,2,j)
%     plot(cell2mat(RMSSD_AS(pat(i),:)),'x b')
%     title(['Patient ' num2str(pat(1,i)) ' with Window size: ' num2str(win(1,j))])
%    hold on
%     plot(cell2mat(RMSSD_QS(pat(i),:)),'x r')
%     hold off
%     ylabel('RMSSD');
%     xlabel('index')
%     legend('AS', 'QS')
%     ylim([0 0.05])
% 
% 
%     end
%      title([ ' with Window size: ' num2str(win(1,j))])
% end