function [loc,z]=outspikeloc(x)
%outputs outspiked points and locations
if size(x,2)>=size(x,1)
    x=x';
end

if (size(x,2)==1)
    loc=(find(x<(mean(x)-2.7*std(x))));
    for i=1:length(loc)
        if loc(i)>11
            x(loc(i))=mean(x(loc(i)-10:loc(i)-1));
        elseif loc(i)<(length(x)-11) 
            
            x(loc(i))=mean(x(loc(i)+1:loc(i)+11));
        end
    end
    
    loc=(find(x>(mean(x)+2.7*std(x))));
    for i=1:length(loc)
        if loc(i)>11
            x(loc(i))=mean(x(loc(i)-10:loc(i)-1));
        elseif loc(i)<(length(x)-11) 
            x(loc(i))=mean(x(loc(i)+1:loc(i)+11));
        end
    end
    z=x;
    
else
    for k=1:size(x,2)
        
        loc=(find(x(:,k)<mean(x(:,k)-2.7*std(x(:,k)))));
        for i=1:length(loc)
            if loc(i)>11
                x(loc(i),k)=mean(x((loc(i)-10):(loc(i)-1),k));
            elseif loc(i)<(length(x)-11) 
                
                x(loc(i),k)=mean(x((loc(i)+1):(loc(i)+11),k)); mean(x(loc(i):loc(i)+10),k);
            end
        end
        
        
        
        
        
        loc=(find(x(:,k)>(mean(x(:,k))+2.7*std(x(:,k)))));
        for i=1:length(loc)
            if loc(i)>11
                x(loc(i),k)=mean(x((loc(i)-10):(loc(i)-1),k));
            elseif loc(i)<(length(x)-11) 
                x(loc(i),k)=mean(x((loc(i)+1):(loc(i)+11),k)); mean(x(loc(i):loc(i)+10),k);
                
            end
        end
        
        z=x;
    end
    
end

    
    
        
        
        
        
        
        
        
        
        
        
        
        
        
    end
    

