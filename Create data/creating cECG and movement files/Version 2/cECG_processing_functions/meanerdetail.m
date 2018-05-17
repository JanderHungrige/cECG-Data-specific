%Option to find average range in a signal
%First start by excluding zeros, assuming rawData is an 8x40000 matrix.

function mranvec=meanerdetail(rawData, chunker)
%chunker is how many chunks you want it divided into
sumdata=sum(rawData);
rawData=rawData(:,find(sumdata));
lendata=size(rawData,2);

chunks=[1:chunker]*floor(lendata/chunker);
chunks=[1 chunks];
for i=1:chunker
ranvec(i,:)=max(rawData(:,chunks(i):chunks(i+1)),[],2)-min(rawData(:,chunks(i):chunks(i+1)),[],2);
end

mranvec=mean(ranvec);

    
