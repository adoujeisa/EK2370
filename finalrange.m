clear all;
close all;
clc;


% Read the audiofile
[y,Fs] = audioread('fmcw_double_opposite.wav');%single_target_range.wav'); 
numRidges = 2; %Number of targets

data = -y(:,1);
sync = -y(:,2); 
% Parameters
c = 3e8;                % Speed of light [m/s]
f_start = 2.408e9;            % Start Frequency [Hz]
f_stop = 2.495e9;             % Stop Frequency  [Hz]
bandwidth = f_stop - f_start; % [Hz]
f_center = (f_start + f_stop)/2;
Tp = 20e-3;                   % Pulse width (up-chirp) [s]
N = Tp * Fs;                  % Number of samples per up-chrip pulse
dr = c / (2*bandwidth);       % Range resolution [m]
max_range = (N * dr)/2;       % Maximum range [m]

%% Define search limits and frequency ranges for range and velocity extraction
fridgeLimitRange = 25;
fridgeLengthRange = round(8 * bandwidth * fridgeLimitRange / c); % Factor 4*2 combined for compactness
freqsRange = linspace(fridgeLimitRange, 0, fridgeLengthRange);

fridgeLimitVel = 10;
fridgeLengthVel = round(8 * f_center * fridgeLimitVel / c); % Factor 4*2 combined for compactness
freqsVel = linspace(fridgeLimitVel, 0, fridgeLengthVel);

% Define range vector
range = linspace(0, max_range, 4 * N);

%% Parse sync data and store indices of positive sync values
sync_pulse = [(sync > 0), (1:length(sync)).']; % Create a matrix with sync square waveform and indices

% Find indices where the sync signal is positive
time_temp = find(sync_pulse(:, 1) == 1);


%% Build the matrix where timesteps are row-wise
up_data_parsed = zeros([], N); % Pre-allocate 
down_data_parsed = zeros([],N);
time = zeros(1,[]);            % Same
k = 1;
for i = 2:(size(sync_pulse)-2*N) 
    if sync_pulse(i,1) == 1 && sync_pulse(i-1) == 0 % First value of a row = first up-chirp value
        up_data_parsed(k,:) = data(i:i+N-1)';
        down_data_parsed(k,:) = data(i+N:i+2*N-1)';
        time(1,k) = sync_pulse(i,2) / Fs;
        k = k + 1;
    end
end

%% MS Clutter Rejection
up_data_parsed = bsxfun(@minus, up_data_parsed, mean(up_data_parsed, 1)); % Subtract column mean to each column
down_data_parsed = bsxfun(@minus, down_data_parsed, mean(down_data_parsed, 1)); % Subtract column mean to each column

%% 3-step MTI
MTI3_up = zeros(size(up_data_parsed));
for t = 3:size(up_data_parsed,1)
    MTI3_up(t,:) = up_data_parsed(t,:) - 2*up_data_parsed(t-1,:) + up_data_parsed(t-2,:);
end

MTI3_down = zeros(size(down_data_parsed));
for t = 3:size(down_data_parsed,1)
    MTI3_down(t,:) = down_data_parsed(t,:) - 2*down_data_parsed(t-1,:) + down_data_parsed(t-2,:);
end

%% FFT Define a function for computing and normalizing the FFT of the input data
computeFFT = @(data) max(20*log10(abs(fft(data, 4 * N, 2))), -1000000);

% Compute and normalize the FFT for up-chirp and down-chirp data
sfft_up = computeFFT(MTI3_up);
sfft_up = sfft_up(:, 1:end/2) - max(sfft_up(:)); % Retain first half and normalize

sfft_down = computeFFT(MTI3_down);
sfft_down = sfft_down(:, 1:end/2) - max(sfft_down(:)); % Retain first half and normalize

% Compute the unnormalized FFT for ridge extraction
sfftFridgeUp = fft(MTI3_up, 4 * N, 2);
sfftFridgeDown = fft(MTI3_down, 4 * N, 2);

%% find fridges
% Define a function to extract ridges using tfridge
extractRidges = @(sfftFridge, rangeLength, freqs, numRidges) ...
    tfridge(rot90(sfftFridge(:, 1:rangeLength)), freqs, 1, 'NumRidges', numRidges);

% Extract ridges for up-chirp and down-chirp data
fridge_up = extractRidges(sfftFridgeUp, fridgeLengthRange, freqsRange, numRidges);
fridge_down = extractRidges(sfftFridgeDown, fridgeLengthRange, freqsRange, numRidges);

% Compute beat frequency, range, Doppler frequency, and velocity
f_beat = (fridge_down + fridge_up) / 2;
rangeExtracted = 15 * (f_beat * c * Tp) / bandwidth;

f_doppler = (fridge_down - fridge_up) / 2;
velExtracted = 20 * (c * f_doppler) / (2 * f_center);



%% Plotting
% Plot the captured data
figure;
plot(data);
title('Captured Data');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

% Plot the Range-Spectrogram using imagesc
figure;
imagesc(range, time, sfft_up, [-35 0]); % Set color limits directly
xlim([0, 25]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
title('Range-Spectrogram');

% Plot the Range vs. Time
figure;
plot(time, rangeExtracted);
title('Range vs Time');
xlabel('Time [s]');
ylabel('Range [m]');
grid on;
