clear all;
close all;
clc;

%% COMPUTE VELOCITY FROM AUDIO FILE
% Read the audio file
[y, fs] = audioread('SAR_Test_File.m4a'); 

% read and store in matrix, and inversion by the sound card
backscatter = -y(:, 1); % Intensity of the received signal
sync = -y(:, 2);

% Parameters
c = 3e8;                % Speed of light [m/s]
f_start = 2.408e9;      % Start Frequency [Hz]
f_stop = 2.495e9;       % Stop Frequency [Hz]
bandwidth = f_stop - f_start; % Bandwidth [Hz]
delta_r = c / (2*bandwidth); % Range resolution [m]
f_center = (f_start + f_stop) / 2; % Center Frequency [Hz]
Tp = 0.02;              % Pulse width [s]
N = round(Tp * fs);     % Number of samples per pulse
N_total = length(backscatter); % Total number of samples in read file
max_range = (N * delta_r)/2; % Maximum range [m]
num_pulse = floor((N_total - N) / (N * 2)); % Total number of pulses

Trp = 0.25;         % Duration of measuring at each position [s]
Nrp = Trp*Fs;       % Number of samples at each position
Tp = 20e-3;         % Upchirp length [s]
N = Tp*Fs;          % Number of samples per upchirp

lambda = c/fc;      % WaveLength [m]
delta_x = lambda/2;
