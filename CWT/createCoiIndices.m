 function indices = createCoiIndices(N)
    if rem(N,2)  % is odd
        indices = 1:ceil(N/2);
        indices = [indices, fliplr(indices(1:end-1))];
    elseif ~rem(N,2)  % is even
        indices = 1:N/2;
        indices = [indices, fliplr(indices)];
    end
 end
