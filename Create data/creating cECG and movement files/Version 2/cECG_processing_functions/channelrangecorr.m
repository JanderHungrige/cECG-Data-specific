%Find which channels have better range
function channelssorted=channelrangecorr(s1,s1p,win1,win2,howmany)
%s1 is the raw signal. s1p is the processed version after Hilbert
%transforming it.
%find the ones with the biggest range, then drop out the badly correlated
%ones and use the top ones (suggest 3) for reconstruction. 
for i=1:8
    ran(i)=range(s1{i}(win1:win2));
end
[rt,et]=sort(ran,'descend');
channels=et(1:howmany);
bestchannel=channels(1);
sens=cell2mat(s1p);
cors=corrcoef(sens);
[val,badchann]=find(cors(bestchannel,:)<0.5);
%  [y,my]=sort(sum(cors));
%  badchann=my(1:2);
count=1;
for i=1:howmany
    if  ~ismember(channels(i),badchann)
        channelssorted(count)=channels(i);
        count=count+1;
    end
end
    
   

