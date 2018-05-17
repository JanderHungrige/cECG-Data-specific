% frequency domain measures of HRV
% this function calculates the frequency domain HRV features over the
% entire spectrum, only divided into the differnt states AS,QS,... while the " frequencydomainHRV.m" calculate it over 30
% sec epoch

%%%%long term analysis (best 24h)
% total power variance of all NN intervals  (<0.4 Hz)
% ULF power in ultra low frequency range (adults <0.003 Hz)
% VLF power in very low frequency range (adult 0.003-0.04 Hz)
% LF  power in low frequency range (adult 0.04-0.15 Hz)
% HF power in high frequency range (adult 0.15-0.4 Hz)
% alpha SLope linear interpolation of the spectrum in a log-log scale (<0.04 Hz)


function freqdomainHRV_single (pat,plotting,win,overlap,saving)

for i=1:length(pat)
        for j=1:length(win)
            for k=1:length(overlap)
                
 %%%%%%%%%%%%%%%%%%%%%%long term analysis (best 24h)

%Loading power spectrum and frequency for AS and QS

                cd('C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\ECG no care separated into AS and QS (total)\HRV analysis\Ralphs\')
                if exist(fullfile(cd, ['PowerspectrumAS_' num2str(pat(1,i)) '.mat']), 'file') == 2 % ==> 0 or 2
                   load(fullfile(cd, ['PowerspectrumAS_' num2str(pat(1,i)) '.mat']))
                end
                if exist(fullfile(cd, ['PowerspectrumQS_' num2str(pat(1,i)) '.mat']), 'file') == 2 % ==> 0 or 2
                   load(fullfile(cd, ['PowerspectrumQS_' num2str(pat(1,i)) '.mat']))
                end
                
                % if files are missing pxxAS and fAS might be reused,therfore we put it in a if clause
                if exist(fullfile(cd, ['PowerspectrumAS_' num2str(pat(1,i)) '.mat']), 'file') == 2 ...
                    && exist(fullfile(cd, ['PowerspectrumQS_' num2str(pat(1,i)) '.mat']), 'file') == 2

    % total power variance of all NN intervals  (<0.4 Hz)

                    totpowAS_tot=sum(pxxAS(fAS<0.4));
                    totpowQS_tot=sum(pxxQS(fQS<0.4));

                    if saving
                        folder='C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\ECG no care separated into AS and QS (total)\HRV analysis\Ralphs\Freqdomain\';
                        save([folder 'totpowAS_tot' num2str(pat(i)) '.mat'],'totpowAS_tot');
                        save([folder 'totpowQS_tot' num2str(pat(i)) '.mat'],'totpowQS_tot');
                    end

    % ULF power in ultra low frequency range (adults <0.003 Hz)

                    ULFAS_tot=sum(pxxAS(fAS<0.003));
                    ULFQS_tot=sum(pxxQS(fQS<0.003));

                    if saving
                        folder='C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\ECG no care separated into AS and QS (total)\HRV analysis\Ralphs\Freqdomain\';
                        save([folder 'ULFAS_tot' num2str(pat(i)) '.mat'],'ULFAS_tot');
                        save([folder 'ULFQS_tot' num2str(pat(i)) '.mat'],'ULFQS_tot');
                    end

    % VLF power in very low frequency range (adult 0.003-0.04 Hz)

                    VLFAS_tot=sum(pxxAS(fAS>0.003 & fAS<0.04));
                    VLFQS_tot=sum(pxxQS(fQS>0.003 & fQS<0.04));

                    if saving
                        folder='C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\ECG no care separated into AS and QS (total)\HRV analysis\Ralphs\Freqdomain\';
                        save([folder 'VLFAS_tot' num2str(pat(i)) '.mat'],'VLFAS_tot');
                        save([folder 'VLFQS_tot' num2str(pat(i)) '.mat'],'VLFQS_tot');
                    end

    % LF  power in low frequency range (adult 0.04-0.15 Hz)

                    LFAS_tot=sum(pxxAS(fAS>0.04 & fAS<0.15));
                    LFQS_tot=sum(pxxQS(fQS>0.04 & fQS<0.15));

                    if saving
                        folder='C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\ECG no care separated into AS and QS (total)\HRV analysis\Ralphs\Freqdomain\';
                        save([folder 'LFAS_tot' num2str(pat(i)) '.mat'],'LFAS_tot');
                        save([folder 'LFQS_tot' num2str(pat(i)) '.mat'],'LFQS_tot');
                    end

    % HF power in high frequency range (adult 0.15-0.4 Hz)

                    HFAS_tot=sum(pxxAS(fAS>0.15 & fAS<0.4));
                    HFQS_tot=sum(pxxQS(fQS>0.15 & fQS<0.4));

                    if saving
                        folder='C:\Users\310122653\Documents\PhD\InnerSense Data\Matlab\saves\ECG no care separated into AS and QS (total)\HRV analysis\Ralphs\Freqdomain\';
                        save([folder 'HFAS_tot' num2str(pat(i)) '.mat'],'HFAS_tot');
                        save([folder 'HFQS_tot' num2str(pat(i)) '.mat'],'HFQS_tot');
                    end

    % alpha Slope linear interpolation of the spectrum in a log-log scale (<0.04 Hz)

                end%if exist
                 clearvars -except pat win i j k saving overlap plotting

            end %overlap
        end%win
end%pat
end
