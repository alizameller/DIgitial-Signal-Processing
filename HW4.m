% DSP HW4
% Aliza Meller
close all; 
clc;
clear;

%% Question 1: Sig Gen
% Part a
fs0 = 2000; 
fs = 16*fs0;
B = 800;
M = 16;

N = 2^(nextpow2(fs));  % Will use N as number of samples to generate
Nfft = 2*N;

binfreqs = colon(0,(fs/Nfft),fs);
[c0, k0] = min(abs(binfreqs - B)); %find index that has a value closest to 800Hz
% Because MATLAB indexes from 1 and not 0, although it is calculated as
% 820, the actual index is at 819
k0 = k0 - 1; % k0 = 819
xf = randn(1, k0 + 1); %% Generates iid N(0,1) samples, which will occupy DFT indices 0<=k<=k0
x = real(ifft(xf, Nfft)); %% Indices k0+1<=k<=Nfft-1 are all 0
x = x(1:N);
x = x/std(x);  % Normalize to unit variance

% Part b
% While xf is white because it has zero mean and is uncorrelated, it is not
% stationary because the variance is not constant over the signal because
% of the zero-padding. x is white and stationary because the inverse 
% fourier transform a white signal (xf) white and stationary. 

% Part c
[pxx, w] = pwelch(x);
f = w*(fs/(2*pi)); 
plot(f,pxx);
xlabel("Frequency (Hz)")
ylabel("Power/Frequency (dB/(rad/sample))")
title("Power Spectral Density of x")
hold on
plot(800, pxx(800), 'r*');

% Part d
% The properties of the fourier transform tell us that the inverse FT of a
% real signal is symmetric and complex. Because xf is real, the ifft
% produces a symmetric and complex signal as the output. 

% Part e 
% Assume ergodicity. It is called unit variance because the variance of
% x/std(x) = 1. See below. 
variance = var(x); 

% Part f
% In order for the samples to be uncorrelated, the PSD needs to be flat. In 
% part c it is clear that the PSD is not flat which neans they must be 
% correlated.  

%% Question 2: Quantizer
samples = randi([-5, 5], 1, 20);
Qbits = 3; 
scale = 3; 
quantized = quantop(samples, Qbits, scale);

% By observation, the upper bound is 2.25 and the lower bound is -3. This is
% because the range went from -3 to 3 with 2^3 (= 8) steps. 6/8 = 0.75, so
% each step is 0.75. 

%% Question 3: Decimation Filter
% Part a
wcrit = pi/M; 
% passband: 0 <= w <= 7/8 wc
% stopband: wc <= w <= pi
passband_edge = (0.875*wcrit*fs)/(2*pi); % convert passband edge into hertz
% passband_edge = 875 > B = 800

% Part b
fcrit = [0, (7/8)*wcrit, wcrit, 1];
% nominal gain in passband: 1, stopband: 0
acrit = [1 1 0 0]; 
[b0, a0] = firpm(221, fcrit, acrit); % returns desired amplitude of freq. 
% response of the resultant filter, b and the maximum ripple height ERR, a
b0 = b0/norm(b0);
a0 = 1; % normalize filter to unit power

% Part c
[H,W] = freqz(b0, a0);
figure; 
plot((W*fs)/(2*pi), abs(H));
title('Magnitude Response of Decimation Filter');
xlabel('Frequency (Hz)');
ylabel('Magnitude (Linear)');

% Part d
% There is amplitude distortion as seen by the rippling in the amplitude 
% and there is no phase distortion because it is a linear FIR filter. 

%% Question 4: Differential Coder and Accumulator
% Part a
% Transfer Function of Accumulator:
% b0 = 0, b1 = 1, a1 = 1 
% H(Z) = (b0 + b1 * z^-1) / (1 - a1 * z^-1)
% H(Z) = z^-1/(1 - z^-1) = 1/(z - 1)

% Transfer Function of Decoder:
% H(Z) = yn/xn = 1/[(1/Acc) + 1] = 1/[(z-1) + 1] = 1/z

% Part b
y = zeros(1, N);
e = zeros(1, N);
e(1) = x(1);

for i = 2:N
    y(i) = e(i-1) + y(i-1);
    e(i) = x(i) - y(i);
end

var_e = var(e);

% Looking at the first 1-20 values of x are the same as the first 2-21 values
% of y which shows a shift by 1

% Part c
% See PDF document 

%% Question 5: Running the Systems
% Part a
% To find alpha, we convert the given statement to P((|X-mu|/sigma) > alpha) = 0.01
% Because X ~ N(mu, sigma^2), Z = |X-mu|/sigma ~ N(0, 1). 
% Using the Z score, we can rewrite the statement as P(|Z| > alpha) < 0.01
% which can be written as P(|Z| > alpha) = P(Z > alpha) + P(Z < -alpha) 
% So, we want to find P(Z > alpha) < 0.005 and P(Z < -alpha) < 0.005
% The function norminv provides the inverse cdf for the normal dist. which
% allows us to find alpha given the paramter 0.995. (or, -alpha given 0.005)
alpha = norminv(0.995);

% Part b
% Unquantized:
unquantized_12 = conv(b0, x);
unquantized_12 = downsample(unquantized_12, M);
% System I:
Qbits1 = 5;
Scale1 = 10;
x1 = quantop(unquantized_12, Qbits1, Scale1);
% System II:
Qbits2 = 5;
Scale2 = 5;
x2 = downsample((conv(quantop(x, Qbits2, Scale2), b0)), M);

figure;
subplot(3,1,1)
pwelch(unquantized_12)
title("Unquantized Spectral Output Sytems 1 and 2")
subplot(3,1,2)
pwelch(x1)
title("Quanitzed Spectral Output System 1")
subplot(3,1,3)
pwelch(x2)
title("Quantized Spectral Output System 2")

% Computing SNR
error1  = unquantized_12 - x1;
error2 = unquantized_12 - x2;
SNR1 = 10*log10(mean(unquantized_12.^2)/mean(error1.^2));
SNR2 = 10*log10(mean(unquantized_12.^2)/mean(error2.^2));

% Part c
% System 3
unquantized_3 = conv(y, b0);
unquantized_3 = downsample(unquantized_3, M);

Qbits3 = 5;
Scale3 = 0.5;

y3 = zeros(1, N);
e3 = zeros(1, N);
e3(1) = x(1); 

for i = 2:N
    y3(i) = quantop(e3(i-1), Qbits3, Scale3) + y3(i-1);
    e3(i) = x(i) - y3(i);
end

x3 = conv(b0, y3);
x3 = downsample(x3, M);
error3 = unquantized_3 - x3;
SNR3 = 10*log10(mean(unquantized_3.^2)/mean(error3.^2));

figure;
nexttile
pwelch(unquantized_3)
title("Unquantized Output of System 3")
nexttile
pwelch(x3)
title("Quantized Output of System 3")

%% Question 6
% Part a
% System III has the greatest SNR, and system I has the smallest SNR.

% All results discussed below correspond to an SNR of system I = 22.40 dB, 
% SNR of system II = 27.51 dB, SNR of system III = 37.78 dB
octaves = log2(M); 
SNRbenefit_per_octave21 = (SNR2 - SNR1)/octaves; % SNR benefit of II over I = 1.28 dB
SNRbenefit_per_octave31 = (SNR3 - SNR1)/octaves; % SNR benefit of III over I = 3.85 dB

% Part b
bits_precision_21 = (SNR2 - SNR1)/6; % 0.85 increased bits of precision of II over I
bits_precision_31 = (SNR3 - SNR1)/6; % 2.56 increased bits of precision of III over I

% Part c
% The SNR benefit of system II over I is 1.28 dB

% Part d
% Our standard assumption of roundoff error is to model it as an additive
% noise. The quantization noise occupies the badwidth of the 
% signal, fs/2. Oversampling by an octave reduces the signal bandwidth to 
% fs/4. Because the PSD of AWGN is flat, the noise still occupies fs/2. 
% As a result, to get rid of the extra noise we filter the signal 
% through a decimation filter to avoid aliasing of higher frequencies. 
% This filter needs to be a low pass filter with a cutoff freq at the dropoff
% point, of the signal which is fs/4 = (fs/2) / 2. This decreases the noise 
% power by a factor of 2, and the overall SNR increases by a factor of 2.  
% Therefore, because 10*log10(2) = 3 dB, SNR increases by 3 dB. 


function xq = quantop(x,Qbits,scale)
    x = x/scale;
    xq = round((x+1)*2^(Qbits-1))/2^(Qbits-1)-1;
    xq = max(min(xq,1-2^-(Qbits-1)),-1); % -1 <= x < 1
    xq = xq*scale;
end



