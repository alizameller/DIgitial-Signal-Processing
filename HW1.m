%Aliza Meller DSP HW1
close all; 
clc;
clear;

%% Question 2
% Part c
fs = 50000; 
binfreqs = colon(0,(fs/256),fs);
f1 = 10000; 
f2 = 40000; 

% index1-1 and index2-1 are the indices where peaks occur 
[c1, index1] = min(abs(binfreqs-f1)); %find index that has a value closest to 10kHz
[c2, index2] = min(abs(binfreqs-f2)); %find index that has a value closest to 40kHz

k1 = binfreqs(index1); 
k2 = binfreqs(index2);

offset1 = abs(k1 - f1); 
offset2 = abs(k2 - f2);

%calculating the straddle loss
straddle1 = 20*log10(1/abs(diric((2*pi/fs)*(offset1), 250))); % not 256 because it is padded with 0s
straddle2 = 20*log10(1/abs(diric((2*pi/fs)*(offset2), 250)));

% Part d
f = 20000;
n = 0:249; 

% obtain the hamming window
hammingWindow = (hamming(250)).';

offset3 = sum(exp(-1*(2*pi*offset1/fs)*j*n) .* hammingWindow); 
offset4 = sum(exp(-1*(2*pi*offset2/fs)*j*n) .* hammingWindow); 

w0 = abs(sum(exp(0*(n)) .* hammingWindow));

%calculate straddle loss
straddle3 = -20*log10((abs(offset3))/w0);
straddle4 = -20*log10((abs(offset4))/w0);

% Part e
binfreqs2 = colon(0,(fs/1024),fs);

[c3, index3] = min(abs(binfreqs2-f1)); 
[c4, index4] = min(abs(binfreqs2-f2));

k3 = binfreqs2(index3);
k4 = binfreqs2(index4);

offset5a = abs(k3 - f1); 
offset6a = abs(k4 - f2);

%computing straddle loss without a hamming window
straddle5a = 20*log10(1/abs(diric((2*pi/fs)*(offset5a), 250)));
straddle6a = 20*log10(1/abs(diric((2*pi/fs)*(offset6a), 250)));

offset5b = sum(exp(-1*(2*pi*offset5a/fs)*j*n) .* hammingWindow); 
offset6b = sum(exp(-1*(2*pi*offset6a/fs)*j*n) .* hammingWindow); 

%computing straddle loss with a hamming window
straddle5b = -20*log10((abs(offset5b))/w0);
straddle6b = -20*log10((abs(offset6b))/w0);

%% Question 4
% Part b
N = 200;
W0 = zeros(1, N);
Wn = ones(1, N);

%for loop used to make W0[4m] = 1 for integer m
for i = 1:200
    if mod(i-1, 4) == 0 
           W0(i) = 1; 
    end
end

%computing 200 point DFTs of W[n] and W0[n]
DFT1_0 = abs(fft(W0, N));
DFT1_n = abs(fft(Wn, N)); 

figure;
subplot(2,1,1);
stem(DFT1_0);
title('DFT Magnitude Spectra for W0[n]');
xlabel('Frequency (Hz)');
ylabel('Magnitude (linear)');

subplot(2,1,2);
stem(DFT1_n);
title('DFT Magnitude Spectra for W[n]');
xlabel('Frequency (Hz)');
ylabel('Magnitude (linear)');

hold on;

% Part c
N0 = 2^12;
%calculating DFTs of rectangular window W[n] and W0[n]
DFT2_0 = fftshift(fft(W0, N0));
DFT2_n = fftshift(fft(Wn, N0)); 

normalizedDFT2_0 = abs(DFT2_0)./max(DFT2_0);
normalizedDFT2_n = abs(DFT2_n)./max(DFT2_n); 

%range in radians
rad = 2*pi*(-(N0)/2:((N0)/2)-1)/(N0);

figure;
subplot(2,1,1);
plot(rad, normalizedDFT2_0);
title('4096 point DFT Magnitude Spectra for W0');
xlabel('Frequency (rad)');
ylabel('Magnitude (linear)');

subplot(2,1,2);
plot(rad, normalizedDFT2_n);
title('4096 point DFT Magnitude Spectra for Wn');
xlabel('Frequency (rad)');
ylabel('Magnitude (linear)');

% Part d
A = 1;
theta = 0;
w_0 = 0.4;
xn = A*cos((w_0 + theta)*(0:199));
xn_W0 = xn.*W0;

%computing DFT of both the rectangular window and W0[n]
DFT3_xn = abs(fft(xn, N0)/N0);
DFT3_xn_W0 = abs(fft(xn_W0, N0)/N0);

%range in radians
rad2 = 2*pi*(0:N0-1)/N0;

figure
plot(rad2, DFT3_xn);
title('X[n] DFT Magnitude Spectra for rectangular window and W0[n]');
hold on
plot(rad2, DFT3_xn_W0);
legend('DFT for rectangular window','DFT for W0[n]')
%setting the range to 0 <= w <= pi
xlim([0, pi]);
xlabel('Frequency (rad)');
ylabel('Magnitude (linear)');

% Part e
% The spectral resolution can be seen through the fuzziness of the plots. In
% this plot, the spectral resolution appears to be the same for the two
% windows. 

%% Question 5
load handel.mat
% Part a
y8s = y(1:65536);
yBlocks = reshape(y8s, Fs, 8);

DFT = abs(fft(yBlocks, Fs));

figure
plot(DFT);
title('DFT Magnitude Spectra (linear scale) of music sample');
xlabel('Frequency (Hz)');
ylabel('Magnitude (linear)');
xlim([0, Fs/2]);
 
% There are peaks at around 574Hz and 1125 Hz and 1125 is approximmately
% twice 574. 
% The frequencies we are hearing are those peak frequences (i.e. the loudest ones) and the 
% harmonics (the smaller peaks) contribute to how we hear the sound even though they are hard
% to make out because they are used to determine the pitch. 

% Part b
pxx = mag2db(pwelch(y, Fs/2, Fs/4, Fs));
figure
plot(pxx);
title('Periodogram for the whole music signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
%limiting the x-range from 0 to Fs/2
xlim([0, Fs/2]);

% Part c
% The peaks in the DFTs of the individual music segments are at the same 
% frequencies as in the periodogram. These peaks occur at approximate integer multiples
% of ~574Hz (~1125Hz, ~1696Hz, ~2260Hz) and they represent the harmonics. 

% Part d
figure
spectrogram(y,Fs/2,Fs/8,Fs/2);
title('A graph of the spectrogram');
xlabel('Frequency (Hz)');
ylabel('Time (s)');


