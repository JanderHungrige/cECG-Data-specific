function outs= flip_proximity(yu,proximity_matrix)
%take a ranking matrix from worst to bset channel and flip accordimng to
%proximity
fu=flipud(yu);
if fu(1)~=8
prox_data=proximity_matrix(fu(1),:);
ranking=[8,7,6,5,4,3,2,1];
for i=1:length(fu)
    weight(i)=ranking(i)*prox_data(fu(i));
end
[iu,iis]=sort(weight,'descend');

outs=[fu(1) fu(iis)'];
outs=outs(1:8);
else 
    outs=fu';
end


