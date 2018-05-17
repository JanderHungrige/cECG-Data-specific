function create_cECG_movement() 
%SCRIPT for Jan Werth
%
% This script extracts the capacitve ECG + R peak positions + a measure of
% the motion level in the mat files called:
% >> 'motion_level.mat' (250 Hz)
% >> 'cap_ECG.mat' (250 Hz) 
% >> 'cap_Rpeakposition.mat' (in number of samples)
% >> 'ref_Rpeakposition.mat' (in number of samples)
%
% NOTE: ALL the signals above are delayed by 1.74 seconds compared to the
% reference ECG !! (also the ref_Rpeakpositions!!!)
%
%  >> 'referenceECG'(250 Hz)
%
%
% aline.serteyn@gmail.com
% 13-4-2017

            launch_processing
            
    %% Display results for one folder
    FS = 500; 
    % figure; aax(1) = subplot(4,1,1); 
    % plot([1:length(cap_ECG)]./FS, cap_ECG); hold on; 
    % plot(cap_Rpeakposition./FS, cap_ECG(cap_Rpeakposition), 'ok'); title('Capacitive ECG and detected peaks')
    % aax(2) = subplot(4,1,2);
    % plot(RRpositions(1:end-1)./FS, diff(RRpositions)./FS,'k'); hold on
    % plot(RRpositionsREF(1:end-1)./FS, diff(RRpositionsREF)./FS,'r'); title('RR tacogram'); ylabel('RR interval [s]'); legend('cap ECG', 'ref ECG'); ylim([0.2 0.9]) %HR 66 tot 300 bpm
    % aax(3) = subplot(4,1,3);
    % plot([1:length(motion_level)]./FS, motion_level); title('Motion level (low motion = below 1e8, high motion = above 1e9)')
    % xlabel('Time [s]')
    % aax(4) = subplot(4,1,4); 
    % ref_RpeakpositionSYNC = max(1,ref_Rpeakposition-FS*1.74);
    % %%% If you have the reference ECG:
    % plot([1:length(referenceECG)]./FS, referenceECG); hold on
    % plot(ref_RpeakpositionSYNC./FS, referenceECG(ref_RpeakpositionSYNC), 'or'); title('Reference ECG and detected peaks')
    % %%% Otherwise:
    % % plot(ref_RpeakpositionSYNC./FS, ones(1,length(ref_RpeakpositionSYNC)), 'or'); title('Reference ECG and detected peaks')
    % xlabel('Time [s]')
    % linkaxes(aax, 'x')
end
toc

