%%        DEPARTMENT OF ELECTRICAL AND COMPUTER ENGINEERING, 
%%                  UNIVERSITY OF THESSALY
%%
%%           CS692: SPEECH AND AUDIO PROCESSING
%%
%% INSTRUCTOR: GERASIMOS POTAMIANOS (gpotamianos@inf.uth.gr)
%%
%% PROJECT BY: NONAS EVANGELOS (vagnonas@gmail.com),
%%             CHATZIGEORGIOU CHRYSOSTOMOS (hrhatzig@gmail.com)
%%

clear; clc;

%% PART B: SIMPLE SYSTEM OF SPEECH PRODUCTION
%%

% Parameters setup.
Fs = 10e3; Ts = 1/Fs;            % sampling frequency and period
total_duration = 2;              % total duration time (sec)
phonem_duration = 0.5;           % phonem duration time (sec)

s = 30;                          % wide bandwidth 2*sk = 60Hz
A = 5;                           % gain

% Average formant frequencies for the vowels of American English.
F = [570, 840,  2410;            %  /AO/
     270, 2290, 3010;            %  /IY/
     300, 870,  2240;            %  /UH/
     530, 1840, 2480];           %  /EH/

pitchPeriod = [8e-3, 4e-3];      % pitch period

%% B6-B10 Create speech signal for given pitch period
%%
for k = 1:length(pitchPeriod)    % fro each pitch period 
    
    Tp = pitchPeriod(k);
    Np = Tp*Fs;
        
    for l = 1 : 4                % for each vowel 
        
        f = F(l,:);              % get formant frequencies for each vowel

        %%  B1) Voiced excitation p[n]

        p_num = (1);
        p_den = zeros(1, Np+1);
        p_den(1) = 1;
        p_den(81) = -0.9999;

        if (k == 1 && l == 1)
            fvtool(p_num, p_den);
        end
        
        [p_sig, ~] = impz(p_num, p_den);
        
        %%  B2) Glotal pulse g[n]
        
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


        %% B3) Vocal tract impulse response v[n]

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
        
        %% B4) Radiation load r[n]

        r = zeros(1,10);
        r(1) = 1;
        r(2) = -0.96;
        
        if (k == 1 && l == 1)
            fvtool(r);
        end
        
        %% B5) Create vowel signal s[n] via convolution s[n] = A(p[n]*g[n]*v[n]*r[n])
        
        h1 = conv(p_sig, g_sig);
        h2 = conv(h1, v_sig);
        h3 = conv(h2, r);
        
        s_sig = A*h3;

       %% Append vowels to create the desired speech signal
       
        if l == 1 
            speech_signal = s_sig;
        else
            speech_signal = [speech_signal; s_sig];
        end
        
    end
    
    % write speech signal to output file
    filename = char(strcat('pitch_', int2str(Np), '.wav'));
    audiowrite(filename, speech_signal, Fs);
    
    % plot speech signal
    figure;
    plot(speech_signal);
    figure_title = sprintf('Speech signal for Np = %d', Np);
    title(figure_title);
end
