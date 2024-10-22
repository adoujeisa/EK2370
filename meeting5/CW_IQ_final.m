close all;
clear all;
clc;


%% COMPUTE VELOCITY FROM AUDIO FILE
% Read the audiofile 
[I,Fs] = audioread('CW_single_re.wav'); 
[Q,Fs] = audioread('CW_single_im.wav'); 
nTargets = 2;

% Composite the complex value and its conjugate
data1 = complex(I,Q);
data2 = conj(data1);

% Parameters
c = 3e8;                % Speed of light [m/s]
f_center = 2.45e9;             % Center Frequency [Hz]
Tp = 0.1;                     % Pulse width [s]
N = Tp * Fs;                  % Number of samples per pulse

% Parse the data
X = mod(-mod(length(data1), N), N);      % Used to find the previous divisible value with respect to length(up_data)
data_cut = data1((N-X+1):end);           % Remove the first elements so that we can reshape up_data
data_parsed = reshape(data_cut, N, [])';
final_data1 = bsxfun(@minus, data_parsed, mean(data_parsed, 2)); % Subtract column mean to each column

X = mod(-mod(length(data2), N), N);      % Used to find the previous divisible value with respect to length(up_data)
data_cut = data2((N-X+1):end);           % Remove the first elements so that we can reshape up_data
data_parsed = reshape(data_cut, N, [])';
final_data2 = bsxfun(@minus, data_parsed, mean(data_parsed, 2)); % MS Clutter rejection

% FFT 
% Compute the FFT, convert to dB scale, and take the first half of the spectrum
f1 = 20 * log10(abs(fft(final_data1, 4 * N, 2)));
f1 = f1(:, 1:end/2);

f2 = 20 * log10(abs(fft(final_data2, 4 * N, 2)));
f2 = f2(:, 1:end/2);

% Define frequency and velocity axes
delta_f = linspace(0, Fs/2, size(f1, 2)); 
vel = (delta_f * c) / (2 * f_center);
time = linspace(1, Tp * size(f1, 1), size(f1, 1));

% Normalize by subtracting the row-wise maximum values
f1_norm = f1 - max(f1, [], 2);
f2_norm = f2 - max(f2, [], 2);



% Find ridges and compute velocities
fridge1 = tfridge(rot90(f1(:, 1:82)), delta_f(1:82), 1, 'NumRidges', nTargets);
fridge2 = tfridge(rot90(f2(:, 1:82)), delta_f(1:82), 1, 'NumRidges', nTargets);

vel1 = (c * fridge1) / (2 * f_center);
vel2 = (c * fridge2) / (2 * f_center);



% Plot
% Plot the first figure
figure;
imagesc(vel, time, f1_norm);
clim([-10 0]);
colorbar;
xlim([0 10]);
xlabel('Velocity [m/sec]');
ylabel('Time [sec]');
title('I + jQ');

% Plot the second figure
figure;
imagesc(vel, time, f2_norm);
clim([-10 0]);
colorbar;
xlim([0 10]);
xlabel('Velocity [m/sec]');
ylabel('Time [sec]');
title('I - jQ');




