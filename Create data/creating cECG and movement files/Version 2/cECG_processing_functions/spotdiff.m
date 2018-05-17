function [S,ord]= spotdiff(filteredData,filteredDataREF)

  y=(filteredDataREF);
  y=(y-mean(y))/std(y);
count=1;  
for i=1:7
    for j=i+1:8
        if i~=j
            x=(filteredData(i,:)-filteredData(j,:));
            x=(x-mean(x))/std(x);
             R=corrcoef(abs(x),abs(y));
            S(count)=abs(R(2,1));
            M(count,:)=[i,j];
            %figure;plot(x);hold all;plot(y); pause;
            
            count=count+1;
            
        end
    end
end
%Now find the maximal correlation values
[val, loc]=max(S);

ord=M(loc,:);

