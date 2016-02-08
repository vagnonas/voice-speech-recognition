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

%% PART A: BASIC ALGORITHMS FOR VOCAL PARAMETERS ESTIMATION

figure_counter = 1;   % counter for figure ploting
RecordAudio = 1;      % change to 0 to read an existing wav file.

%% A1) Create a recorder object to record speaker's name.
%%

disp('A1) Recording...');

Fs = 16e3;      % sample rate
bps = 8;        % bits per sample
ca = 1;         % chanel audio (mono)
rt = 3;         % record time is 3 seconds

if RecordAudio == 1
    recorder = audiorecorder(Fs,bps,ca);

    disp('Start speaking.')
    recordblocking(recorder, rt);
    disp('End of Recording.');

    x = getaudiodata(recorder);
    audiowrite('./myname.wav',x,Fs);
else
    disp('Getting audio data from existing file...');
    [x, ~] = audioread('./myname.wav');
end
disp('Press any key to continue...');
pause;

%% A2) Spectrogram using Hamming short-time and longer-time window 
%%

disp('A2) Spectrogram using Hamming short-time and longer-time window');

fft_size = 4096;

wl = 10e-3;           % window time duration
wo = 5e-3;            % window time overlap

L = floor(wl * Fs);   % window length
overlap = (wo/wl)*L;  % window overlap

% plot spectrogram of short-time window
figure(figure_counter);
spectrogram(x, hamming(L), L-overlap, fft_size, 'yaxis');
title('Spectrogram using short-time Hamming window');

figure_counter = figure_counter + 1;

wl = 100e-3;          % window time duration
wo = 5e-3;            % window time overlap

L = floor(wl * Fs);   % window length
overlap = (wo/wl)*L;  % window overlap

% plot spectrogram of longer-time window
figure(figure_counter);
spectrogram(x,hamming(L),overlap);
title('Spectrogram using longer-time Hamming window');

figure_counter = figure_counter + 1;

disp('Press any key to continue...');
pause;

%% A3) Voiced/unvoived speech segments estimation
%%

disp('A3) Voiced/unvoived speech segments estimation');

wl = 30e-3;          % window time duration
wo = 15e-3;          % window time overlap

L = floor(wl * Fs);  % window length
overlap = (wo/wl)*L; % window overlap

w = hamming(L);      % define hamming window of length L

Eth = 1.5;           % define short-time energy threshold
Zth = 0.1;           % define zero-crossing rate threshold

step = (L-overlap);
stop = (length(x)-L);

for k = 1 : step : stop
    xshort = x(k:k+L-1);
    
    [En, Zn] = extract_features(xshort, w);
    
    energy(k:k+L-1) = En;
    zcr(k:k+L-1) = Zn;
    
    if (En > Eth) && (Zn < Zth)     % speech segment is voiced
        segment(k:k+L-1) = 2;
    elseif (En < Eth) && (Zn > Zth) % speech segment is unvoiced
        segment(k:k+L-1) = 1;
    else                            % speech segment is silence
        segment(k:k+L-1) = 0;
    end
end

% plot sugnal, short-time energy & zero-crossing rate
figure(figure_counter);
subplot(3,1,1);
plot(x);
xlabel('n');
ylabel('x[n]');
title('Speech signal');

figure(figure_counter);
subplot(3,1,2);
plot(energy);
xlabel('n');
ylabel('En');
title('Short-time energy');

figure(figure_counter);
subplot(3,1,3);
plot(zcr);
xlabel('n');
ylabel('Zn');
title('Zero-crossing rate');

figure_counter = figure_counter + 1;

% plot recorded signal & decision levels
figure(figure_counter);
subplot(2,1,1);
plot(x);
xlabel('n');
ylabel('x[n]');
title('Speech signal');

figure(figure_counter);
subplot(2,1,2);
plot(segment);
xlabel('n');
y_values = [0,1,2];
y_labels = {'silence', 'unvoiced', 'voiced'};
set(gca, 'Ytick',y_values,'YTickLabel',y_labels);
title('Decision levels');

figure_counter = figure_counter + 1;

disp('Press any key to continue...');
pause;

%% A4) Pitch period estimation with autocorrelation values & cepstrum
%%

disp('A4) Pitch period estimation with autocorrelation values & cepstrum');

% Pitch period estimation with autocorrelation signal.

R = zeros(1, length(segment));   % auto-correlation signal
ppr = zeros(1, length(segment)); % pitch period in each segment

prit_segment = 1;
for k = 1 : step : stop
    xshort = x(k:k+L-1);
    
    if segment(k) == 2   % if segment is voiced
        r = autocorrelation(xshort, w);
        pp = pitch_period_estimation(r);
        
        % plot autocorrelation of the first voiced segment
        if prit_segment == 1
            figure(figure_counter);
            plot(r);
            title('Autocorrelation signal for a specific segment');
            figure_counter = figure_counter + 1;
            prit_segment = 0;
        end
        
        R(k:k+L-1) = r;
        ppr(k:k+L-1) = pp;
    end
end

% plot autocorrelation signal & pitch period
figure(figure_counter);
subplot(3,1,1);
plot(x);
title('Speech signal');

figure(figure_counter);
subplot(3,1,2);
plot(R);
title('Autocorrelation signal of voiced segments');

figure(figure_counter);
subplot(3,1,3);
plot(ppr);
title('Pitch period estimated with auto-correlation method');
figure_counter = figure_counter + 1;

% Pitch period estimation using cepstrum.

cep = zeros(1, length(segment)); % cepstrum signal
ppc = zeros(1, length(segment)); % pitch period in each segment

prit_segment = 1;
offset = 50; len = 200;
for k = 1 : step : stop
    xshort = x(k:k+L-1);
    
    if segment(k) == 2   % if segment is voiced
        c = cepstrum(xshort);
        cep(k:k+L-1) = c;
        
        % plot cepstrum of the first voiced segment
        if prit_segment == 1
            figure(figure_counter);
            plot(c);
            title('Cepstrum of a specific voiced segment');
            figure_counter = figure_counter + 1;
            prit_segment = 0;
        end
        
        [~, ipos] = max(c(offset:len));
        pp = offset+ipos;      
        ppc(k:k+L-1) = pp;
    end
end

% plot cepstrum signal & pitch period

figure(figure_counter);
subplot(3,1,1);
plot(x);
title('Speech signal');

figure(figure_counter);
subplot(3,1,2);
plot(cep);
title('Cepstrum of voiced segments');

figure(figure_counter);
subplot(3,1,3);
plot(ppc);
title('Pitch period estimated with cepstrum method');

figure_counter = figure_counter + 1;

% plot speech signals, pitch period estimated with autocorrelation and
% cepstrum

figure(figure_counter);
subplot(3,1,1);
plot(x);
title('Speech signal');

figure(figure_counter);
subplot(3,1,2);
plot(ppr);
title('Pitch period of voiced segment using autocorrelation');

figure(figure_counter);
subplot(3,1,3);
plot(ppc);
title('Pitch period of voiced segments using cepstrum');

figure_counter = figure_counter + 1;

disp('Press any key to continue...');
pause;

%% A5) Linear Prediction
%%

disp('A5) Linear prediction');

S = zeros(2, L); % S(1, :) is the unvoiced segment,
                 % S(2, :) is the voiced segmrnt

% find first unvoiced segment
for k = 1 : step : stop    
    if segment(k) == 1
        S(1,:) = x(k:k+L-1);
        break;
    end
end

% find first voiced segment
for k = 1 : step : stop    
    if segment(k) == 2
        S(2,:)= x(k:k+L-1);
        break;
    end
end

% linear prediction
for p = 8 : 4 : 16         % for each filter order
    fprintf(' order p = %d\n', p);
    for k = 1:2            % for voiced/unvoiced segments
        sn = S(k,:);
    
        rn = autocorrelation(sn, w);     % segment autocorrelation
        [M, rhsvec] = toeplitz(rn, p);   % toplitz matrix & rhs vector
        a = M\rhsvec;                    % system solution
        e = lpc_error(sn, w, a, p);      % prediction error
        
        if k == 1 % for unvoiced segment
            fprintf('   unvoiced segment prediction error = %f\n', e);
        else      % for voiced segment
            fprintf('   voiced segment prediction error = %f\n', e);
        end
        
        % plot all-pole trasfer function magnitude and roots on Z-plane
        [b,a] = zp2tf([],a,1);
        fvtool(b,a);
        
        % plot signal segment dft
        figure(figure_counter);
        plot(abs(fft(sn)));
        if k == 1 % for unvoiced segment
            title('Unvoiced segment DFT');
        else      % for voiced segment
            title('Voiced segment DFT');
        end
        figure_counter = figure_counter + 1;
        
        if p == 8   % plot voiced/unvoiced segment dft the first time only
            % plot signal segment dft
            figure(figure_counter);
            plot(abs(fft(sn)));
            if k == 1 % for unvoiced segment
                title('Unvoiced segment DFT');
            else      % for voiced segment
                title('Voiced segment DFT');
            end
            figure_counter = figure_counter + 1;
        end
    end
end