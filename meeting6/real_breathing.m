clear all;
close all;
clc;

% Read the real and imaginary parts of the signal
[I, Fs] = audioread('SDR_CWIF_BREATHING_JONNE2_REAL.wav'); 
%[Q, ~] = audioread('SDR_CWIF_BREATHING_JONNE2_IMAG.wav'); 

% Given constants
c = 3e8; % Speed of light (m/s)
f0 = 2.4e9; % Example frequency (Hz) for a 2.4 GHz signal


data = I(100000:end, 1);
% Calculate the phase angle phi between I and Q
%phi = atan2(Q, I); % atan2 provides the four-quadrant inverse tangent

% Calculate the phase difference over time (to handle wrapping)
%dphi = diff(unwrap(phi)); % Unwrap the phase to avoid discontinuities

% Calculate the respiratory rate using the formula
Ri = c * data / (4 * pi * f0);

% Generate a time vector adjusted for the size of Ri
t_Ri = (0:length(Ri)-1) / Fs;

% Plot the respiratory rate over time
figure;
plot(t_Ri, Ri);
xlabel('Time (s)');
ylabel('Respiratory Rate (m/s)');
xlim([0, 20]); % 设置 x 轴范围

title('Respiratory Rate over Time');
grid on;

