function e = lpc_error(s, w, a, p)
    e = 0;
    for m = p+1 : length(s)-p-1
        sum = 0;
        for k = 1 : p
            sum =  sum + a(k)*s(m-k);
        end
        e = e + (s(m)-sum)^2;
    end
end

