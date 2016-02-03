%% DEPARTMENT OF ELECTRICAL AND COMPUTER ENGINEERING, 
%%             UNIVERSITY OF THESSALY
%%
%%        HY692: SPEECH AND AUDIO PROCESSING
%%
%% INSTRUCTOR: GERASIMOS POTAMIANOS (gpotamianos@inf.uth.gr)
%% PROJECT BY: NONAS EVANGELOS (vagnonas@gmail.com),
%%             CHATZIGEORGIOU CHRYSOSTOMOS (hrhatzig@gmail.com)
%%

clear; clc;

%% PART B: SIMPLE SYSTEM OF SPEECH PRODUCTION

% parameters setup;
Fs = 10e3; Ts = 1/Fs;    % sampling frequency and period
total_duration = 2;      % total duration time (sec)
phonem_duration = 0.5;   % phonem duration time (sec)

s = 30;                  % wide bandwidth 2*sk = 60Hz
A = 5000;                % gain

% Average formant frequencies for the vowels of American English.
F = [570, 840,  2410;    %  /AO/
     270, 2290, 3010;    %  /IY/
     300, 870,  2240;    %  /UH/
     530, 1840, 2480];   %  /EH/

pitchPeriod = [8e-3, 4e-3];      % pitch period

%% B1) 
for k = 1:length(pitchPeriod)
    
    Tp = pitchPeriod(k);
    Np = Tp*Fs;
    
    for l = 1 : 4
        
        f = F(l,:);              % formant frequencies for each vowel

        %%  Voiced excitation p[n]
%         n = (0:1023);
%         p = zeros(1, 1024);
%         for i = 1 : 80 : length(p)
%             p(i) = 0.9999^i;
%         end
%         figure(1);
%         plot(n,p);
%         title('Voiced Excitation (time-domain)');
%         xlabel('n');
%         ylabel('p[n]');

        p_num = (1);
        p_den = zeros(1, 81);
        p_den(1) = 1;
        p_den(81) = -0.9999;

        if (k == 1 && l == 1)
            fvtool(p_num, p_den);
        end
        
        [p_sig, Np] = impz(p_num, p_den);
        
        %% Glotal pulse g[n]
        
        Ng = (0:39);
        g_sig = zeros(1, 40);
        for i = 1 : 25
            g_sig(i) = 0.5*(1-cos(pi*(Ng(i)+1)/25));
        end
        for i = 26 : 34
            g_sig(i) = cos(0.5*pi*(Ng(i)-24)/10);
        end

        if (k == 1 && l == 1)
            fvtool(g_sig);
        end


        %% Vocal tract impulse response v[n]

        v_num = (1);
        v_den = (1);
        for i = 1 : 3
          denpoly = [1, -2*exp(-2*pi*s*Ts)*cos(2*pi*f(i)*Ts), exp(-4*pi*s*Ts)];
          v_den = conv(v_den, denpoly);
        end


        if (k == 1 && l == 1)
            fvtool(v_num, v_den);
        end
        
        [v_sig, Nv] = impz(v_num, v_den);
        
        %% Radiation load r[n]

        r = zeros(1,10);
        r(1) = 1;
        r(2) = -0.96;
        
        if (k == 1 && l == 1)
            fvtool(r);
        end
        
        %% Create vowel signal s[n] via convolution s[n] = A(p[n]*g[n]*v[n]*r[n])
        
        h1 = conv(p_sig, g_sig);
        h2 = conv(h1, v_sig);
        h3 = conv(h2, r);
        
        s_sig = A.*h3;
        plot(s_sig);
        
        break;
    end
    break;
end
