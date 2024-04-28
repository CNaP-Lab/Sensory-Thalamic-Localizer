function [out] = removeSlice(in, direction, numToShrink)
    
    out = in;
    if ( strcmpi(direction,'superior') )
        plane = false(size(in,1),size(in,2),1);
        dimension = 3;
        for i = 1:numToShrink
            anyOne = squeeze(any(out,[1,2]));
            firstOne = find(anyOne,1,'last');
            out(:,:,firstOne) = plane;
        end

    elseif ( strcmpi(direction,'posterior') )
        plane = false(size(in,1),1,size(in,3));
        dimension = 2;
        for i = 1:numToShrink
            anyOne = squeeze(any(out,[1,3]));
            firstOne = find(anyOne,1,'first');
            out(:,firstOne,:) = plane;
        end

    elseif ( strcmpi(direction,'inferior') )
        plane = false(size(in,1),size(in,2),1);
        dimension = 3;
        for i = 1:numToShrink
            anyOne = squeeze(any(out,[1,2]));
            firstOne = find(anyOne,1,'first');
            out(:,:,firstOne) = plane;
        end
        
    else
        pause(eps); drawnow;
        error('Invalid direction specified');
    end
    
end
