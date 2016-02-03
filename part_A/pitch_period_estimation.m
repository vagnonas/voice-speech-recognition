function [pp] = pitch_period_estimation(R)
    
    [~, i] = max(R);
    i1 = i;
    [~, j] = min(R(i1:length(R)));
    j1 = j;
    
    [~, i] = max(R(j1:length(R)));
    i2 = i + j1;
    
    [~, j] = min(R(i2:length(R)));
    j2 = j + i2;
    
    [~, i] = max(R(j2:length(R)));
    i3 = i + j2;
        
    pp = 0.5*(i2-i1)+0.5*(i3-i1);

end

