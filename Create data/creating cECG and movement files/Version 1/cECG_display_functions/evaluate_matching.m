
function saved_length=evaluate_matching(fid, error, tolerance, WINDOW_SIZE)

% Get the % of time when there is a good match with REF
timeCov=round(sum(error<=tolerance)/length(error)*100); % NOTE: 1 meanRR per processing window (4000) => 2 Hz!

% Count number of unreliable periods longer than 10s
lengthy=0; count=0;
saved_length=[];
for i=1:length(error)
    if error(i)>tolerance
        lengthy=lengthy+1; %in number of processing windows
    else
        if lengthy>10*8000/WINDOW_SIZE %>10 seconds
            count = count+1;
            saved_length=[saved_length, lengthy]; 
        end
        lengthy=0;
    end
end
if lengthy>10*8000/WINDOW_SIZE %>10 seconds
    count = count+1;
    saved_length=[saved_length, lengthy]; 
end
saved_length=saved_length./8000*WINDOW_SIZE; % in seconds!

if ~isempty(fid)
    fprintf(fid,'%s\t\t\t%2.0f\t\t\t%2.0f\t\t\t%2.1f\n',[num2str(tolerance) ' bpm'], timeCov, count, max(saved_length));
end