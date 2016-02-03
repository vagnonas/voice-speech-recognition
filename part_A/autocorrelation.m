function [r] = autocorrelation(x, w)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = length(x);
r = zeros(1, N);

for k = 1 : N
    for m = 1 : N-k
        r(k) = r(k) + (x(m)*w(m)) * (x(m+k)*w(m+k)); 
    end
end


end

