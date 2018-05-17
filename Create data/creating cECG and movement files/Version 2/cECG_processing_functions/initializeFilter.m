
function filterCoefficients=initializeFilter(Fs,cutoff,order,type)


if strcmp(type,'lowpass')==1
    [b,a]=butter(order,2*cutoff/Fs,'low');
elseif strcmp(type,'highpass')==1
    [b,a]=butter(order,2*cutoff/Fs,'high');
elseif strcmp(type,'bandstop')==1
    if length(cutoff)==1
        [b,a]=butter(order,2*cutoff/Fs,'stop');
    else
        error('provide 2 values for bandstop cutoff');
    end
else
    error('unknown filter type');
end
   

%Calculate filter specific variables
b=b(:); a=a(:);
na = numel(a);
nb = numel(b);
L = 1;
% Check coefficients

nfilt = max(nb,na);
nfact = 3*(nfilt-1);  % length of edge transients
% Zero pad shorter coefficient vector as needed
if nb < nfilt
    b(nfilt,1)=0;
elseif na < nfilt
    a(nfilt,1)=0;
end

% Compute initial conditions to remove DC offset at beginning and end of
% filtered sequence.  Use sparse matrix to solve linear system for initial
% conditions zi, which is the vector of states for the filter b(z)/a(z) in
% the state-space formulation of the filter.
if nfilt>1
    rows = [1:nfilt-1, 2:nfilt-1, 1:nfilt-2];
    cols = [ones(1,nfilt-1), 2:nfilt-1, 2:nfilt-1];
    vals = [1+a(2), a(3:nfilt).', ones(1,nfilt-2), -ones(1,nfilt-2)];
    rhs  = b(2:nfilt) - b(1)*a(2:nfilt);
    zi   = sparse(rows,cols,vals) \ rhs;   
else
    zi = zeros(0,1);
end

filterCoefficients=struct('b',b,'a',a,'zi',zi,'nfact',nfact);