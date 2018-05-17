% This function determines the indices of the R & S peaks in the ECG and
% the  index at which the RS slope is steepest.
%
% Inputs:
% t_ecg = ECG time axis (s)
% ecg   = ECG [V]
% fs    = Sampling frequency [Hz]
%
% Outputs:
% ecg_r_peak_idx   = array holding the R peak indices
% ecg_s_peak_idx   = array holding the S peak indices
% ecg_qr_slope_idx = array holding the QS slope indices
% ecg_rs_slope_idx = array holding the RS slope indices
% ecg_st_slope_idx = array holding the ST slope indices
% bbi_ecg          = array of beat-to-beat intervals in ms
% bpm_ecg          = array of BPM

%by Wijshoff, Ralph <ralph.wijshoff@philips.com>

function  Resp_find_peaks(t_Resp, Resp, fs, Resp_max, ploting,saving)

% Determine the derivative of the respiration
diff_resp = fs*[0; diff(Resp)*(-1)];
threshold = 1*mean(diff_resp(diff_resp > 0));

% Determine the threshold to select the sharply decreasing RS slope
rs_threshold = 1*mean(diff_resp(diff_resp < 0));

% Find the indices of all RS slopes
rs_slope_idxs = find(diff_resp < rs_threshold);

% Select the minimum of each of the RS slopes,
% so the index of the steepest decrease of the RS slope
rs_slope_start = [1; 1+find(diff(rs_slope_idxs)>1)];
num_beats      = length(rs_slope_start);
Resp_qr_slope_idx = zeros(1, num_beats);
Resp_rs_slope_idx = zeros(1, num_beats);
Resp_st_slope_idx = zeros(1, num_beats);
Resp_r_peak_idx   = zeros(1, num_beats);
Resp_s_peak_idx   = zeros(1, num_beats);
Nqrs           = fs*(50e-3);
Nrec           = length(Resp);

for i = 1:num_beats
    % Select the interval in which the RS slope is to be found
    if i < num_beats
        intv = rs_slope_idxs(rs_slope_start(i)) : ...
               rs_slope_idxs(rs_slope_start(i+1) - 1);
    else
        intv = rs_slope_idxs(rs_slope_start(i):end);
    end
    
    % Determine the index at which the RS slope is steepest
    [tmp, min_idx]      = min(diff_resp(intv));
    Resp_rs_slope_idx(i) = min_idx + intv(1) - 1;
    
    % Determine the index at which the QR slope is steepest
    intv = (Resp_rs_slope_idx(i) - Nqrs) : Resp_rs_slope_idx(i)+fs*0.25;
    intv(intv<1) = [];
    [tmp, max_idx]      = max(diff_resp(intv));
    Resp_qr_slope_idx(i) = max_idx + intv(1) - 1;
    
    % Determine the index at which the ST slope is steepest
    intv = Resp_rs_slope_idx(i) : (Resp_rs_slope_idx(i) + Nqrs+fs*0.25);
    intv(intv>Nrec) = [];
    [tmp, max_idx]      = max(diff_resp(intv));
    Resp_st_slope_idx(i) = max_idx + intv(1) - 1;    
    
    % Determine the maximum of the ECG in this interval 
    % to identify the R peak
    intv = Resp_qr_slope_idx(i) : Resp_rs_slope_idx(i)+fs*0.25;
    [tmp, max_idx] = max(Resp(intv));
    Resp_r_peak_idx(i) = max_idx + intv(1) - 1;
    
    % Determine the minimum of the ECG in this interval 
    % to identify the S peak
    intv = Resp_rs_slope_idx(i) : Resp_st_slope_idx(i)+fs*0.25;
    [tmp, min_idx] = min(Resp(intv));
    Resp_s_peak_idx(i) = min_idx + intv(1) - 1;    
end

% Remove RS slopes and the corresponding R&S peaks that deviate too much
% from the average.
% Test2 2011/11/10 @ 4 km/h: subtract 0.1 from the deviations
qr_slope_deviation     = -3.5;
rs_slope_deviation     = 1.1;
st_slope_deviation     = 2.5;
rs_amplitude_deviation = 1.4;
rs_amplitudes          = Resp(Resp_r_peak_idx) - Resp(Resp_s_peak_idx);

qr_slope_limit = qr_slope_deviation*mean(diff_resp(Resp_qr_slope_idx));
rs_slope_limit = rs_slope_deviation*mean(diff_resp(Resp_rs_slope_idx));
st_slope_limit = st_slope_deviation*mean(diff_resp(Resp_st_slope_idx));
rs_peak_limit  = rs_amplitude_deviation*mean(rs_amplitudes);

not_ok         = (diff_resp(Resp_rs_slope_idx) > rs_slope_limit) | ...
                 ( (diff_resp(Resp_qr_slope_idx) < qr_slope_limit) & ...
                   (diff_resp(Resp_st_slope_idx) < st_slope_limit) ) | ...
                 (rs_amplitudes < rs_peak_limit);

Resp_qr_slope_idx(not_ok) = [];
Resp_rs_slope_idx(not_ok) = [];
Resp_st_slope_idx(not_ok) = [];
Resp_r_peak_idx(not_ok)   = [];
Resp_s_peak_idx(not_ok)   = [];

% Also remove any peaks which have been found twice
% (which may be a result of motion artifacts)
not_ok = find(diff(Resp_r_peak_idx)==0);
Resp_qr_slope_idx(not_ok) = [];
Resp_rs_slope_idx(not_ok) = [];
Resp_st_slope_idx(not_ok) = [];
Resp_r_peak_idx(not_ok)   = [];
Resp_s_peak_idx(not_ok)   = [];

% Remove any spurious peaks
nmin = round((60/Resp_max)*fs);
idx  = find(diff(Resp_r_peak_idx) < nmin, 1);
while(~isempty(idx))
    % Remove the one having the smallest R peak
    % [tmp del_idx] = min(ecg(ecg_r_peak_idx([idx idx+1])));
    % Remove the one having the more shallow slope
    [tmp del_idx] = max(diff_resp(Resp_rs_slope_idx([idx idx+1])));
    
    not_ok = idx + (del_idx - 1);
    
    Resp_qr_slope_idx(not_ok) = [];
    Resp_rs_slope_idx(not_ok) = [];
    Resp_st_slope_idx(not_ok) = [];
    Resp_r_peak_idx(not_ok)   = [];
    Resp_s_peak_idx(not_ok)   = [];
    
    idx  = find(diff(Resp_r_peak_idx) < nmin, 1);
end

% Correct 

% Determine the beat-to-beat-intervals [ms] and BPM
bbi_Resp = [nan, 1000*(diff(Resp_r_peak_idx))/fs];
bpm_Resp = 60e3 ./ bbi_Resp;

    if ploting
        figure('Name', 'Resp peak detection'), 
        h1 = subplot(3,1,1); hold all, 
            plot(t_Resp, Resp, 'LineWidth', 1),
            plot(t_Resp(Resp_r_peak_idx), Resp(Resp_r_peak_idx), '*r'),
            plot(t_Resp(Resp_s_peak_idx), Resp(Resp_s_peak_idx), '*c'),
            plot(t_Resp(Resp_qr_slope_idx), Resp(Resp_qr_slope_idx), '*g'),
            plot(t_Resp(Resp_rs_slope_idx), Resp(Resp_rs_slope_idx), '*m'),
            plot(t_Resp(Resp_st_slope_idx), Resp(Resp_st_slope_idx), '*y'),
            grid on, box on, axis tight,
            xlabel('Time [s]', 'FontSize', 12),
            ylabel('Amplitude [V]', 'FontSize', 12),
            title([' Resp; number of beats: ', num2str(num_beats)], 'FontSize', 12),
            legend('Resp', 'R peak', 'S peak', 'QR slope', 'RS slope', 'ST slope')
            set(gca, 'FontSize', 12),
        h2 = subplot(3,1,2); hold all,
            plot(t_Resp, diff_resp, 'LineWidth', 1),
            plot([t_Resp(1) t_Resp(end)], rs_threshold*[1 1], 'm'),
            plot([t_Resp(1) t_Resp(end)], rs_slope_limit*[1 1], '--m'),
            plot([t_Resp(1) t_Resp(end)], qr_slope_limit*[1 1], '--g'),
            plot([t_Resp(1) t_Resp(end)], st_slope_limit*[1 1], '--y'),
            plot(t_Resp(Resp_qr_slope_idx), diff_resp(Resp_qr_slope_idx), '*g'),
            plot(t_Resp(Resp_rs_slope_idx), diff_resp(Resp_rs_slope_idx), '*m'),
            plot(t_Resp(Resp_st_slope_idx), diff_resp(Resp_st_slope_idx), '*y'),
            grid on, box on, axis tight,
            xlabel('Time [s]', 'FontSize', 12),
            ylabel('Amplitude [V/s]', 'FontSize', 12),
            title(' Resp: derivate w.r.t. time', 'FontSize', 12),
            legend('d/dt Resp', 'RS slope threshold', 'RS slope limit', ...
            'QR slope limit', 'ST slope limit', 'QR slope', 'RS slope', 'ST slope')
            set(gca, 'FontSize', 12),
        h3 = subplot(3,1,3); hold all, 
            plot(t_Resp(Resp_r_peak_idx), bpm_Resp, '.b', 'MarkerSize', 10),
            grid on, box on, axis tight,
            xlabel('Time [s]', 'FontSize', 12),
            ylabel('Breathing Rate [BPM]', 'FontSize', 12),
            title('Peak derived respiration', 'FontSize', 12),
            set(gca, 'FontSize', 12),
        linkaxes([h1 h2 h3], 'x');
    end


return