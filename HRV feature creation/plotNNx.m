function plotNNx(plotting, win, pat, overlap)
if plotting
    for k=1:length(pat)
        clearvars NN50_AS NN30_AS NN20_AS NN10_AS NN50_QS NN30_QS NN20_QS NN10_QS

        for j=1%:length(win)
             for o=1%:length(overlap)
               % clear NN50_AS NN30_AS NN20_AS NN10_AS NN50_QS NN30_QS NN20_QS NN10_QS
                 cd('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\windowed_ECG_chunks\HRV analysis\Ralphs\')
                 if exist(fullfile(cd, ['NNx' num2str(pat(1,k)) 'win' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,o)) '.mat']), 'file') == 2 % ==> 0 or 2
                   load((fullfile(cd, ['NNx' num2str(pat(1,k)) 'win' num2str(win(1,j)) 's_win_overlap_' num2str(overlap(1,o)) '.mat'])))
                   
                   %getting rifd of zeroes for plotting
                   NN50_AS(NN50_AS==0)=nan;
                   NN30_AS(NN30_AS==0)=nan;
                   NN20_AS(NN20_AS==0)=nan;
                   NN10_AS(NN10_AS==0)=nan;
                   
                   NN50_QS(NN50_QS==0)=nan;
                   NN30_QS(NN30_QS==0)=nan;
                   NN20_QS(NN20_QS==0)=nan;
                   NN10_QS(NN10_QS==0)=nan;

                   %if in ms bring it to s
                   NN50_AS./1000;
                   NN30_AS./1000;
                   NN20_AS./1000;
                   NN10_AS./1000;
                   
                   NN50_QS./1000;
                   NN30_QS./1000;
                   NN20_QS./1000;
                   NN10_QS./1000;
                                 
                  
                figure(1)
                   plot(NN50_AS,'x b')
                   hold on 
                   plot(NN50_QS,'x r')
                   %ylim([0 0.5])
                  % x=[3 4 6 7 9 13];
                 %  set(gca,'XTick',x); % Change x-axis ticks
                   title(['NN50'])
                   ylabel('NN50')
                   xlabel('index')
               figure(2)
                   plot(NN30_AS,'x b')
                   hold on 
                   plot(NN30_QS,'x r')
                   %ylim([0 0.5])
                   %x=[3 4 6 7 9 13];
                   %set(gca,'XTick',x); % Change x-axis ticks
                   title(['NN30'])
                   ylabel('NN30')
                   xlabel('index')
               figure(3)
                   plot(NN20_AS,'x b')
                   hold on 
                   plot(NN20_QS,'x r')
                   
                   %ylim([0 0.5])
                   %x=[3 4 6 7 9 13];
                   %set(gca,'XTick',x); % Change x-axis ticks
                   title(['NN20'])
                   ylabel('NN20')
                   xlabel('index')
               figure(4)
                   plot(NN10_AS,'x b')
                   hold on 
                   plot(NN10_QS,'x r')
                   
                   %ylim([0 0.5])
                   %x=[3 4 6 7 9 13];
                   %set(gca,'XTick',x); % Change x-axis ticks
                   title(['NN10'])
                   ylabel('NN10')
                   xlabel('index')
               
                figure(5)
                   plot(pat(1,k),NN50_AS','x b')
                   hold on
                   plot(pat(1,k)+0.3,NN50_QS','x r')
                   title(['NN50 over allpatients'])
                   ylabel('NN50')
                   xlabel('patients')
                   xlim([2 14])
                figure(6)
                   plot(pat(1,k),NN10_AS','x b')
                   hold on
                   plot(pat(1,k)+0.3,NN10_QS','x r')
                   title(['NN10 over allpatients'])
                   ylabel('NN10')
                   xlabel('patients')
                   xlim([2 14])
                   
                figure(7)
                 subplot(221)
                   plot(pat(1,k),NN50_AS','x b')
                   hold on
                   plot(pat(1,k)+0.3,NN50_QS','x r')
                   title(['NN50 over allpatients'])
                   ylabel('NN50')
                   xlabel('patients')
                   xlim([2 14])
                 subplot(222)
                   plot(pat(1,k),NN30_AS','x b')
                   hold on
                   plot(pat(1,k)+0.3,NN30_QS','x r')
                   title(['NN30 over allpatients'])
                   ylabel('NN30')
                   xlabel('patients')
                   xlim([2 14])
                 subplot(223)
                   plot(pat(1,k),NN20_AS','x b')
                   hold on
                   plot(pat(1,k)+0.3,NN20_QS','x r')
                   title(['NN20 over allpatients'])
                   ylabel('NN20')
                   xlabel('patients')
                   xlim([2 14])
                 subplot(224)
                   plot(pat(1,k),NN10_AS','x b')
                   hold on
                   plot(pat(1,k)+0.3,NN10_QS','x r')
                   title(['NN10 over allpatients'])
                   ylabel('NN10')
                   xlabel('patients')
                   xlim([2 14])
                   
               figure(8)
                 subplot(221)
                   plot(NN50_AS, 'x b')
                   hold on
                   plot(NN50_QS, 'x r')
                   title(['NN50 "time" (idx)'])
                   ylabel('NN50')
                   xlabel('index')
                   %xlim([2 14])
                 subplot(222)
                   plot(NN30_AS, 'x b')
                   hold on
                   plot(NN30_QS, 'x r')
                   title(['NN30 "time" (idx)'])
                   ylabel('NN30')
                   xlabel('index')
                  % xlim([2 14])
                 subplot(223)
                   plot(NN20_AS, 'x b')
                   hold on
                   plot(NN20_QS, 'x r')
                   title(['NN20 "time" (idx)'])
                   ylabel('NN20')
                   xlabel('index')
                   %xlim([2 14])
                 subplot(224)
                   plot(NN10_AS,'x b')
                   hold on
                   plot(NN10_QS,'x r')
                   title(['NN10 "time" (idx)'])
                   ylabel('NN10')
                   xlabel('index')
                   %xlim([2 14])
                   

                end
             end
        end
        
     end
   end
end