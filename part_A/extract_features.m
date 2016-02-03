function [En, Zn] = extract_features(x, w)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% 
% length(x)
% length(w)
% pause;


En = 0;
Zn = 0;

En = En + (x(1)^2)*w(1);
for k = 2 : length(x)
    En = En + (x(k)^2)*w(k);
    Zn = Zn + abs(sign(x(k))-sign(x(k-1)))*w(k);
end

En = log(En);
Zn = Zn./(2*length(w));


end

