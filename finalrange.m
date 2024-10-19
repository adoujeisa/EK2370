clear all;
close all;
clc;


% Read the audiofile
[y,Fs] = audioread('fmcw_double_opposite.wav'); 
numRidges = 2; % Amout of targets

% inversion by the sound card
data = -y(:,1); % Radar backscatter data (received reflected signal)
sync = -y(:,2); % Sync data (square waveform)

% Parameters
c = 3e8;                % Speed of light [m/s]
f_start = 2.408e9;            % Start Frequency [Hz]
f_stop = 2.495e9;             % Stop Frequency  [Hz]
bandwidth = f_stop - f_start; % [Hz]
dr = c / (2*bandwidth);       % Range resolution [m]
Tp = 20e-3;                   % Pulse width [s]
N = Tp * Fs;                  % Number of samples per pulse
max_range = (N * dr)/2;       % Maximum range [m]

fridgeLimit = 25;
fridgeLength = round((2*bandwidth*fridgeLimit*4)/c); % To limit the search of tfridge
freqs = linspace(fridgeLimit, 0, fridgeLength);
range = linspace(0, max_range, 4*N);

% Parse up-chirp data according to the sync data
sync_pulse = [(sync > 0), (1:length(sync)).']; % Create a matrix with sync square waveform and indexes

% Find indices where the sync signal is positive
time_temp = find(sync > 0);


% Build the matrix where timesteps are row-wise
up_data_parsed = zeros([], N); % Pre-allocate 
time = zeros(1,[]);            % Same
k = 1;
for i = 2:(size(sync_pulse)-N) 
    if sync_pulse(i,1) == 1 && sync_pulse(i-1) == 0 % First value of a row = first up-chirp value
        up_data_parsed(k,:) = data(i:i+N-1)';
        time(1,k) = sync_pulse(i,2) / Fs;
        k = k + 1;
    end
end

% MS Clutter Rejection
final_data = bsxfun(@minus, up_data_parsed, mean(up_data_parsed, 1)); % Subtract column mean to each column

%% 3-step MTI using vectorized indexing
MTI3 = [zeros(2, size(up_data_parsed, 2)); up_data_parsed(3:end, :) - 2 * up_data_parsed(2:end-1, :) + up_data_parsed(1:end-2, :)];


%% Compute the FFT of MTI3 with zero-padding to 4*N points along the 2nd dimension
sfft_MTI3 = abs(fft(MTI3, 4 * N, 2));

% Convert to decibels and normalize
sfft_MTI3 = 20 * log10(sfft_MTI3);
sfft_MTI3 = sfft_MTI3(:, 1:end/2);  % Retain only the first half (up to Nyquist)

% Normalize the data to a maximum value of 0 dB
sfft_MTI3 = sfft_MTI3 - max(sfft_MTI3(:));

% Thresholding to avoid -Inf values (set a lower limit)
sfft_MTI3 = max(sfft_MTI3, -1000000);

%% Compute the FFT of MTI3 with zero-padding to 4*N points along the 2nd dimension
sfft_fridge = fft(MTI3, 4 * N, 2);

% Transpose and extract ridges using tfridge
sfft_rotated = rot90(sfft_fridge(:, 1:fridgeLength));
[fridge, ~, ~] = tfridge(sfft_rotated, freqs, 1, 'NumRidges', numRidges, 'NumFrequencyBins', 10);

% Convert extracted ridges to range
rangeExtracted = (15 * fridge * c * Tp) / bandwidth;



%% Plot the Range-Spectrogram using imagesc
figure;
imagesc(range, time, sfft_MTI3, [-35 0]); % Directly set color limits in imagesc
xlim([0 25]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
title('Range-Spectrogram');

% Plot the Range vs Time using plot
figure;
plot(time, rangeExtracted);
xlabel('Time [s]');
ylabel('Range [m]');
title('Range vs Time plot');
