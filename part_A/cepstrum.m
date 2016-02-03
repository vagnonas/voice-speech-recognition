function [cepstrum] = cepstrum(x)
    cepstrum = ifft(log(abs(fft(x))));
end

