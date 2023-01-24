% DSP HW3
% Aliza Meller and Haley Dave 
close all; 
clc;
clear;

for i = 1:2 % Put everything in a loop to avoid repeated code
N = 4*i;  % loops for N = 4 and N = 8
wname = ['db', int2str(N)];
[h0,h1,f0,f1] = wfilters(wname);

%% Part a
figure; 
tiledlayout('flow');
ax1 = nexttile; 
stem(h0);
title("Analysis Lowpass Filter h0 (N = " + N + ")");
ax2 = nexttile; 
stem(h1);
title("Analysis Highpass Filter h1 (N = " + N + ")");
ax3 = nexttile; 
stem(f0);
title("Synthesis Lowpass Filter f0 (N = " + N + ")");
ax4 = nexttile; 
stem(f1); 
title("Synthesis Highpass Filter f1 (N = " + N + ")");

% From the stem plots, it is clear that h0 is the paraconjugate of f0 and h1 is
% the paraconjugate of f1 because the plots of the synthesis polyphase matrices and
% the analysis polyphase matrices are reflections of each other. 
% As a result, we know that these two filters satisfy the paraunitary PR
% property of the filter bank. 

%% Part b
%% Question 1
% H1(z) = (z^-7)H̃0(-z)
% F0(z) = (z^-7)H̃0(z) 
% F1(z) = (z^-7)H̃1(z) = H0(-z) 
%% Question 2
[H0, w0] = freqz(h0, 1);
[H1, w1] = freqz(h1, 1);

figure; 
plot(w0, abs(H0));
hold on;
plot(w1, abs(H1));
hold off;
title("Superimposition of |H1(w)| and |H0(w)| (N = " + N + ")");
legend('|H0(w)|', '|H1(w)|'); 
ylabel('|H(w)|');
xlabel('w (rad)');
%% Question 3
maxdiff = (max(abs(H0).^2 + abs(H1).^2)); 
if round(maxdiff - 2, 1) == 0 % rounding to nearest whole number
    verified = true;
else 
    verified = false;
end 
%% Question 4
E00 = h0(1:2:end);
E01 = h0(2:2:end);
E10 = h1(1:2:end);
E11 = h1(2:2:end);

E = cell(2);
E{1, 1} = E00; 
E{1, 2} = E01;
E{2, 1} = E10;
E{2, 2} = E11; 

R00 = f0(2:2:end);
R01 = f1(2:2:end);
R10 = f0(1:2:end); 
R11 = f1(1:2:end);

R = cell(2);
R{1, 1} = R00; 
R{1, 2} = R01;
R{2, 1} = R10;
R{2, 2} = R11; 

%% Question 5
Epara = transpose(E);
for i = 1 : 2
    for j = 1 : 2
    Epara{i, j} = fliplr(Epara{i, j});
    end
end

% Case N = 4: 
%   N0 = N - 1 = 3
%   c = 1 (by inspection, it is clear that there is no scaling factor)
% Case N = 8:
%   N0 = N - 1 = 7
%   c = 1 (by inspection, it is clear that there is no scaling factor)

%% Question 6
REproduct = cell(2);
REproduct{1, 1} = conv(E00, R00) + conv(E01, R10);
REproduct{1, 2} = conv(E00, R01) + conv(E01, R11);
REproduct{2, 1} = conv(E10, R00) + conv(E11, R10);
REproduct{2, 2} = conv(E10, R01) + conv(E11, R11);

%% Question 7
T = 0.5*(conv(f0, h0) + conv(f1, h1)); 
h0negz = zeros(1, 2*N); 
h1negz = zeros(1, 2*N); 
h0negz(1:2:end) = h0(1:2:end);
h0negz(2:2:end) = -h0(2:2:end);
h1negz(1:2:end) = h1(1:2:end);
h1negz(2:2:end) = -h1(2:2:end);

A = 0.5*(conv(f0, h0negz) + conv(f1, h1negz));
% Finding Error
Ideal_T = zeros(1, 4*N - 1);
Ideal_T(2*N) = 1; 
maxT_error = max(abs(Ideal_T - T));

Ideal_A = zeros(1, 4*N - 1);
maxA_error = max(abs(Ideal_A - A)); 

%% Question 8
firstDeriv = diff(abs(H0).^2)/(w0(2) - w0(1));
secondDeriv = diff(diff(abs(H0).^2))/(w0(2) - w0(1))^2; 
% pad with zeroes because diff removes an element
firstDeriv(512, 1) = 0;
secondDeriv(511,1) = 0;
secondDeriv(512,1) = 0;

figure; 
tiledlayout('flow');
ax1 = nexttile;
plot(w0, firstDeriv);
title("First Derivative of H0(w) (N = " + N + ")");
xlabel("w [rad]");
ylabel("Magnitude of H0(w)");

ax2 = nexttile;
plot(w0, secondDeriv);
title("Second Derivative of H0(w) (N = " + N + ")");
xlabel("w [rad]");
ylabel("Magnitude of H0(w)");

%% Question 9
    if N == 4
        g = zeros(4, 23); 
        c = {h0 h0; h0 h1; h1 h0; h1 h1};
        figure; 
        for j = 1:4
            g(j, :) = conv(c{j, 1}, upsample(c{j, 2}, 2));
            [Gj, wj] = freqz(g(j, :), 1);
            plot(wj, abs(Gj));
            hold on
        end
     
        hold off
        title("N = 4 Magnitude Response of G Filters of 2 Level Tree Structure")
        xlabel("w (rad)")
        ylabel("|G(w)|")
        legend("G0", "G1", "G2", "G3")
    end
    if N == 8
        g = zeros(7, 110); 
        c = {h0 h0 h0; h0 h0 h1; h0 h1 h0; h0 h1 h1; h1 h0 h0; h1 h0 h1; h1 h1 h0; h1 h1 h1};
        figure; 
        for j = 1:8
            g(j, :) = conv(conv(c{j, 1}, upsample(c{j, 2}, 2)), upsample(c{j, 3}, 4));
            [Gj, wj] = freqz(g(j, :), 1);
            plot(wj, abs(Gj));
            hold on
        end
        hold off
        title("N = 8 Magnitude Response of G Filters of 3 Level Tree Structure")
        xlabel("w (rad)")
        ylabel("|G(w)|")
        legend("G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7")
    end
end
