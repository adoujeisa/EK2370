clear all;
close all;
clc;

% Read the audiofile  
[y, Fs] = audioread('SDR_CWIF_BREATHING_JONNE2_REAL.wav'); 

data = y(100000:end, 1);

% Parameters
c = 3e8;                % Speed of light [m/s]
f_center = 2.4e9;             % Center Frequency [Hz]

% Processing
range = data * c / (4 * pi * f_center); 


% Plotting range
figure();
plot(range); 
%xlim([0, 60]); % Set the x-axis limit here

title("Respiratory Rate"); 
grid on;


