
function [previousBlock,dataBlock]=iirFilter(previousBlock, dataBlock, filterCoefficients,iteration,overlap)

for ch=1:size(dataBlock,2)
    %First pass
    if strcmp(iteration,'first')==0
        xt=previousBlock(:,ch);
    else
        xt = -dataBlock(filterCoefficients.nfact+1:-1:2,ch) + 2*dataBlock(1,ch);
    end
    [temp,zo] = filter(filterCoefficients.b,filterCoefficients.a,xt,filterCoefficients.zi*xt(1)); % yc1 not needed
    [yc2,zo] = filter(filterCoefficients.b,filterCoefficients.a, dataBlock(:,ch), zo);
    xt = -dataBlock(end-1:-1:end-filterCoefficients.nfact,ch) + 2*dataBlock(end,ch);
    yc3 = filter(filterCoefficients.b,filterCoefficients.a, xt, zo);
    
    [temp,zo] = filter(filterCoefficients.b,filterCoefficients.a, yc3(end:-1:1), filterCoefficients.zi*yc3(end));
    yc5 = filter(filterCoefficients.b,filterCoefficients.a, yc2(end:-1:1), zo);
    
    y = yc5(end:-1:1);
    
    %previousBlock(:,ch)=dataBlock(end-2*overlap+1:end-overlap,ch);
    previousBlock(:,ch)=dataBlock(1:end-overlap,ch);
    
    %Combine blocks of filtered data
    if strcmp(iteration,'first')==0
        dataBlock(:,ch)=y; %overwrite last part of data from previous iteration (to avoid problems due to reverse filtering)
    else
        dataBlock(:,ch)=y;
    end
end

