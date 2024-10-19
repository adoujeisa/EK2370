clear all;
close all;
clc;


[y,Fs] = audioread('cw_double_opposite.wav'); 
numRidges = 2; % Number of targets

% inversion by the sound card
data = -y(:,1);


% Parameters
c = 3e8;                % Speed of light [m/s]
f_center = 2.45e9;            % Center Frequency [Hz]
Tp = 0.1;                     % Pulse width [s]
N = Tp * Fs;                  % Number of samples per pulse

fridgeLimit = 10;
fridgeLength = round((2*f_center*fridgeLimit*4)/c); % To limit the search of tfridge
freqs = linspace(fridgeLimit, 0, fridgeLength);


% Parse the data
X = mod(-mod(length(data), N), N);      % Used to find the previous divisible value with respect to length(up_data)
data_cut = data((N-X+1):end);           % Remove the first elements so that we can reshape up_data
data_parsed = reshape(data_cut, N, [])';

% MS Clutter Rejection
final_data = bsxfun(@minus, data_parsed, mean(data_parsed, 2)); % Subtract column mean to each column

% Perform FFT with zero-padding to 4*N points along the 2nd dimension
f = abs(fft(final_data, 4*N, 2));
% Convert the amplitude to decibels (dB)
f = 20 * log10(f);
% Retain only the first half of the FFT results
f = f(:, 1:end/2);


% Normalization 1
f1 = f - max(max(f));
delta_f1 = linspace(0, Fs/2, size(f1, 2)); 
vel1 = (delta_f1 * c)/(2 * f_center);
time1 = linspace(1, Tp * size(f1, 1), size(f1, 1));

% Normalization 2
f_row_max = max(f, [], 2);
f2 = f - f_row_max;
delta_f2 = linspace(0, Fs/2, size(f2, 2)); 
vel2 = (delta_f2 * c)/(2 * f_center);
time2 = linspace(1, Tp * size(f2,1), size(f2,1));

% find fridges
sfft_fridge = fft(final_data, 4*N, 2);
[fridge, ~, ~] = tfridge(rot90(sfft_fridge(:, 1:fridgeLength)), freqs, 1,'numRidges',numRidges,'NumFrequencyBins',10);
% Convert to velocity
velExtracted = 160*(fridge * c) / (2*f_center);


%% Plotting
figure();
sgtitle('single target CW velocity measurement');


imagesc(vel1, time1, f1);
clim([-35 0]);
colorbar;
set(gca,'XLim',[0 10]);
xlabel('Velocity [m/sec]'); ylabel('Time [sec]');
title("Velocity spectogram");

figure();
rangePlot = plot(time1,velExtracted);
title("Velocity vs Time plot"); xlabel('Time [s]'); ylabel('Velocity [m/s]');