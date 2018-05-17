function SampEn_QSE(RR,Neonate,saving,savefolder,win,faktor,Session,S)

% #0 clc std 
% #1 calculate matche sample Entropy (SempEn) with range of r from 0.05 to 0.3 *std (20% is in literatur but not always optimal) 
% #2 Do a moving average over each Entropy over r curve to smooth outliers
% #3 Do a moving average over each column of several 5min epoch of the Entropy over r curves. 
    % Thereby, each 5min Entropy over r curve is smoothened and outliers
    % are removed (Low pass Filter)  
% #4 Find the turning point of the Entropy over r curve. 
    %There should be a sudden increase of matches/ decrease of Entropy with rising r which then should strive against a value   
% #5 Find the minimum Entropy   
% #6 Cut out the part between the Minimum and the truning point to create a mean entropy Value value
    % Thereby, we ensure, that we get rid of the beginning extra high Entropy throug too low r 
% #7 find the Entropy that is most similar/close to the meaned Framed entropy 
% #8 find the r of closest Entropy to framed mean Entropy
% #9 use the found optimal r to calculate the optimal SampEn with m=2
% #10 Calculate the Quadratic Sample Entropy

% #0 
for i=1:length(RR)
     SDNN{1,i}=nanstd(RR{i,1});
end

% upsampling  RR is 300s here
% for i=1:length(RR)
%     oldFS=length(RR{i,1})/300;
% %     tmpRR=interp1(RR{i,1},);
% %     RR{i,1}=
% end
%     


% #1
dim=2; 
tau=1;
fact=0.05:0.05:0.3; % the threshold r 5%*sdt to 50%*std in 5% steps
for j=1:length(RR)
    % calculating the matchcount for each r per 5min(30sepoch centered)
    for k=1:length(fact)% going through r
         r=SDNN{1,i}*fact(k);
         
         SampEn_collection(j,k)=SampEntropy(dim, r, RR{j,1}, tau);
    end
end

% #2-#3
% Now use moving average on the r curves of each 5 min to smoothen the curves 
% Calculate the mean per r value toreceive for each 5min epoch a
% smoothened r curve
if length(RR)>20
    winlen=length(RR)/10;
else % If RR is shorter than 20 then just moving average over 2 windows, otherwise over lenght(RR)/10
    winlen=2;
end
SampEn_collection=movmean(SampEn_collection,3,2);  % smoothen each the curve by moving average (mean); 2 means per row 
Meaned_SampEn=movmean(SampEn_collection,winlen,1 ); % smoothen each the curve by moving average of all other curves (actually window of winlen); 1 means per column 

% #4 % From https://nl.mathworks.com/matlabcentral/answers/250257-find-turning-point-in-data
for J=1:size(Meaned_SampEn,1) % Now find Turning point and minimum SempE or maximum matches, Mean between those points
    x = [1 size(Meaned_SampEn(J,:),2)];
    y = [Meaned_SampEn(J,1) Meaned_SampEn(J,end)]; % Two endpoints on the curve "RR";
    % The slope of the line connecting the two endpoints
    m = ( y(2) - y(1) )/( x(2) - x(1) );
    pm= - 1 / m;
    % Point on the curve (xc,yc), point on the line (xl,yl)
    perpDist = zeros(size(Meaned_SampEn,1),1);
    for i = 1:size(Meaned_SampEn,1)
        xc = i ; yc = Meaned_SampEn(i);
        yl = ( (m * xc) + (m^2 * yc) - (m * x(1)) + y(1) )/(1+ m^2);
        xl = xc - m*(yl - yc);
        % distance^2
        d2 = (xl - xc)^2 + (yl - yc)^2;
        perpDist(i) = d2;
    end
  
    [val_turn(J,1), idx_turn(J,1)] = max(perpDist); % get turning point
% #5    
    [val_min(J,1),indx_min(J,1)]=min(Meaned_SampEn(J,:));% get minimal sample entropy
% #6    
    Framed_SampEn{J,1}=Meaned_SampEn(J,idx_turn(J,1):indx_min(J,1)); % Get the frame of SempEn from the turning point to the minimal entropy
    Framed_meaned_SampEn(J,1)=mean(Framed_SampEn{J,1}); %Get the mean SampEn from the framed part
% #7   
    [val_opt_entr(J,1),idx(J,1)] = min(abs(Meaned_SampEn(J,1)-Framed_meaned_SampEn(J,1)));% find the index where the SampEn is most similar to the Framed_meaned_SampEn
end

% #8 
for i=1:length(RR)
    r_opt(i,1)=SDNN{1,i}*fact(idx(i,1)); % Create the optimal r value for the 5min epoch based on the most similar SampEn to the Framed_meaned_SampEn
    r_opt(i,2)=r_opt(i,1)/SDNN{1,i}; % get the % values for comparison
end

for j=1:length(RR)
% #9    
    SampEn(1,j)=SampEntropy(dim, r_opt(j,1), RR{j,1}, tau);
% #10    
    QSE(1,j)=SampEn(1,j)+log(2*r_opt(j,1));
end



if saving
    Saving(SampEn,savefolder,Neonate,win,Session,S) 
    Saving(QSE,savefolder,Neonate,win,Session,S) 
end

end
%% Nested functions

function SampEn=SampEntropy( dim, r, data, tau )
    % SAMPEN Sample Entropy
%   calculates the sample entropy of a given time series data

%   SampEn is conceptually similar to approximate entropy (ApEn), but has
%   following differences:
%       1) SampEn does not count self-matching. The possible trouble of
%       having log(0) is avoided by taking logarithm at the latest step.
%       2) SampEn does not depend on the datasize as much as ApEn does. The
%       comparison is shown in the graph that is uploaded.

%   dim     : embedded dimension
%   r       : tolerance (typically 0.2 * std)
%   data    : time-series data
%   tau     : delay time for downsampling (user can omit this, in which case
%             the default value is 1)
%
%---------------------------------------------------------------------
% coded by Kijoon Lee,  kjlee@ntu.edu.sg
% Mar 21, 2012
% edited by Jan Werth, 2017
%---------------------------------------------------------------------

if nargin < 4, tau = 1; end
if tau > 1, data = downsample(data, tau); end
% tau=1;

for j=1:length(data)
    N = length(data);
    correl = zeros(1,2);
    dataMat = zeros(dim+1,N-dim);
    for i = 1:dim+1
        dataMat(i,:) = data(i:N-dim+i-1);
    end

    for m = dim:dim+1
        count = zeros(1,N-dim);
        tempMat = dataMat(1:m,:);

        for i = 1:N-m
            % calculate Chebyshev distance, excluding self-matching case
            dist = max(abs(tempMat(:,i+1:N-dim) - repmat(tempMat(:,i),1,N-dim-i)));

            % calculate Heaviside function of the distance
            % User can change it to any other function
            % for modified sample entropy (mSampEn) calculation
            D = (dist < r);

            count(i) = sum(D)/(N-dim);
        end
        % collect all MAtches over the different r (k) for each 5min epoch (j)
%         Matches(j,k)=sum(count) % Really, all matches from B and submatches A???? Or oly B or A?   
        correl(m-dim+1) = sum(count)/(N-dim);
    end

    SampEn = log(correl(1)/correl(2)); %correl1=B correl2=A
    
    clearvars dataMat count tempMat correl
end %for lengh data


end

%% Nested saving
    function Saving(Feature,savefolder, Neonate, win,Session,S)
        if exist('Feature','var')==1
            name=inputname(1); % variable name of function input
            save([savefolder name '_Session_' num2str(S) '_win_' num2str(win) '_' Session],'Feature')
        else
            disp(['saving of ' name ' not possible'])
        end       
    end
 