function blockMat = blockify(M, n, idSwitch)
    M2 = M(:);
    M3 = repmat(M2, n,1);
    for j = 1:n
        M3(j:n:end) = M2;
    
    end
    M4 = reshape(M3, n*size(M,1), size(M,2));   
    R = repmat(n,1, size(M,1));
    p = cumsum(R)+1; % where does the next column start
    idx = zeros(1,p(end)); % initialize
    idx(p) = 1; % set the positions to 1
    idx = idx(1:end-1) ; % ignore the last one
    idx = cumsum(idx);
    idx = idx+1;
    blockMat = M4(:,idx);
    if idSwitch
        for i = 1:size(M,1)
            blockMat(((n*(i-1))+1): (n*i),((n*(i-1))+1): (n*i)) = eye(n);
        end
    end
    
end