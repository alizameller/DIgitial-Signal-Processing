%Aliza Meller DSP HW2
close all; 
clc;
clear;

%% Question 1
% Part a
Hznum = [1.0688, 0.4505, 0.4594, 1.0688]; 
Hzden = 8*[1, -1.5055, 1.2630, -0.3778];

z1 = roots(Hznum);
p1 = roots(Hzden);

tiledlayout(2,1)
ax1 = nexttile;
zplane(z1, p1); 
title(ax1,'Pole-Zero Plot of H(z)')

% Part b
Hznum2a = [-0.4954, 1];
Hzden2a = [1, -0.4954];
Hznum2b = [.7626, -1.0101, 1];
Hzden2b = [1, -1.0101, .7626];

b = conv(Hznum2a, Hzden2b) + ... 
    conv(Hzden2a, Hznum2b);
a = 2*conv(Hzden2a, Hzden2b); 

ax2 = nexttile;
zplane(b, a);
title(ax2,'Pole-Zero Plot of H(z) expressed as the sum of two allpass factors')

% The pole zero plots look identitcal which implies that 
% the functions match. 

% Part c
figure(2)
freqz(Hznum, Hzden);
title('Magnitude and Phase Response of H(z)')

% Part d
figure(3)
phasez(Hznum2a, 2*Hzden2a);
hold on; 
phasez(Hznum2b, 2*Hzden2b);
title('Phase Responses of the two all-pass factors of H(z)')
legend('first factor', 'second factor'); 
hold off;

% Part e
HznumRound = round(Hznum * 16)/16;
HzdenRound = round(Hzden * 16)/16;

Hznum2aRound = round(Hznum2a * 16)/16;
Hzden2aRound = round(Hzden2a * 16)/16;
Hznum2bRound = round(Hznum2b * 16)/16;
Hzden2bRound = round(Hzden2b * 16)/16;

z1Round = roots(HznumRound);
p1Round = roots(8*HzdenRound); 

figure(4)
tiledlayout(2,1)
ax3 = nexttile; 
zplane(z1Round, p1Round);
title(ax3,'Pole-Zero Plots of rounded H(z) (to nearest 1/16)')

bRound = conv(Hznum2aRound, Hzden2bRound) + ... 
    conv(Hzden2aRound, Hznum2bRound);
aRound = 2*conv(Hzden2aRound, Hzden2bRound); 

ax4 = nexttile;
zplane(bRound, aRound);
title(ax4,'Pole-Zero Plots of  H(z) rounded to nearest 1/16 expressed as the sum of two allpass factors')

% The rounded coefficients do slightly affect the 
% pole-zero plots, but not enough to move any poles outside 
% the unit circle and therefore do not result in instability. 

% Part f
[H1, w] = freqz(Hznum, Hzden); 
H1dB = 20*log10(abs(H1));
H2 = freqz(HznumRound, HzdenRound);
H2dB = 20*log10(abs(H2)); 
H3 = freqz(bRound, aRound); 
H3dB = 20*log10(abs(H3));

figure(5)
plot(w, H1dB); 
hold on
plot(w, H2dB); 
hold on
plot(w, H3dB);
ylabel('Magnitude (dB)')
xlabel('Normalized Frequency')
title('Magnitude Responses of Systems (Rounding = 1/16)')
legend('Original System','Original System Rounded to 1/16', 'All-Pass System Rounded to 1/16')

% Part g
HznumRound4 = round(Hznum * 4)/4;
HzdenRound4 = round(Hzden * 4)/4;

Hznum2aRound4 = round(Hznum2a * 4)/4;
Hzden2aRound4 = round(Hzden2a * 4)/4;
Hznum2bRound4 = round(Hznum2b * 4)/4;
Hzden2bRound4 = round(Hzden2b * 4)/4;

bRound4 = conv(Hznum2aRound4, Hzden2bRound4) + ... 
    conv(Hzden2aRound4, Hznum2bRound4);
aRound4 = 2*conv(Hzden2aRound4, Hzden2bRound4); 

H14 = freqz(Hznum, Hzden); 
H1dB4 = 20*log10(abs(H1));
H24 = freqz(HznumRound4, HzdenRound4);
H2dB4 = 20*log10(abs(H2)); 
H34 = freqz(bRound4, aRound4); 
H3dB4 = 20*log10(abs(H3));

figure(6)
plot(w, H1dB4); 
hold on
plot(w, H2dB4); 
hold on
plot(w, H3dB4);
ylabel('Magnitude (dB)')
xlabel('Normalized Frequency')
title('Magnitude Responses of Systems (Rounding = 1/4)')
legend('Original System','Original System Rounded to 1/4', 'All-Pass System Rounded to 1/4')

%% Question 2
% Part 1
Fnum = [1, -0.3];
Fden = [0.3, -1];
hf = freqz(Fnum, Fden); 
H = filter_elliptic; 
[num, den] = tf(H); 
Gnum = polyval(num, hf);
Gden = polyval(den, hf);
G = Gnum./Gden;

Hf = freqz(num, den);
HfdB = 20*log10(abs(Hf));

figure(7)
tiledlayout('flow')
ax5 = nexttile;
plot(w, HfdB);
ylabel('Magnitude (dB)')
xlabel('Normalized Frequency')
title('Original (H)')

GdB = 20*log10(abs(G));
ax6 = nexttile;
plot(w, GdB);
ylabel('Magnitude (dB)')
xlabel('Normalized Frequency')
title('Transformed Filter (G)')
% The transformed filter is a high pass filter

wstop = w(GdB <= -30); 
wstop = wstop(end); 
wpass = w(GdB >= -2); 
wpass = wpass(1); 

% Part b
neghf = freqz(-Fnum, Fden);
negGnum = polyval(num, neghf);
negGden = polyval(den, neghf);
negG = negGnum./negGden;

negGdB = 20*log10(abs(negG));
ax7 = nexttile;
plot(w, negGdB);
ylabel('Magnitude (dB)')
xlabel('Normalized Frequency')
title('Transformed Filter (G) using -F(z)')

% Part d
fnum = conv([1, -0.3],[1, 0.4]);
fden = conv([-0.3, 1],[0.4, 1]);

hf2 = freqz(fnum, fden);

Gnum2 = polyval(num, hf2);
Gden2 = polyval(den, hf2);
G2 = Gnum2./Gden2;

GdB2 = 20*log10(abs(G2));
ax8 = nexttile;
plot(w, GdB2);
ylabel('Magnitude (dB)')
xlabel('Normalized Frequency')
title('Transformed Filter (G) using transformed F(z)')

neghf2 = freqz(-fnum, fden);

negGnum2 = polyval(num, neghf2);
negGden2 = polyval(den, neghf2);
negG2 = negGnum2./negGden2;

negGdB2 = 20*log10(abs(negG2));
ax9 = nexttile;
plot(w, negGdB2);
ylabel('Magnitude (dB)')
xlabel('Normalized Frequency')
title('Transformed Filter (G) using transformed -F(z)')

% The transformed filter using transformed F(z) is a band stop filter
% and the transformed filter using transformed -F(z) is a band pass filter
% The order of G is 2 * the order of the prototype, so it is order 4

