function [M, v] = toeplitz(R, p)

    M = zeros(p,p);
    v = zeros(p,1);
    
    for i = 1:p
        M(i,i) = R(1);
        k = 1;
        for j = i+1:p
            M(i,j) = R(k);
            M(j,i) = R(k);
            k = k+1;
        end
        
        v(i) = R(i+1);
    end

end

