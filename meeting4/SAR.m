close all; clear; clc;

% Load the audio data from the file
[y, Fs] = audioread('SAR_Test_File.m4a'); 

%% DATA MATRIX SETUP
% Extract sync and radar backscatter data while correcting for inversion
backscatter = -y(:,1); % Radar backscatter (reflected signal)
sync = -y(:,2); % Sync signal (square waveform)

c = 3e8;                        % Speed of light [m/s]
f_start = 2.408e9;              % Start frequency [Hz]
f_stop = 2.495e9;               % Stop frequency [Hz]
T_rp = 0.25;         % Measurement time for each position [s]
N_rp = T_rp * Fs; % Samples for each position
chirpDuration = 20e-3;            % Length of upchirp [s]
samplesPerChirp = chirpDuration * Fs; % Samples per upchirp
fc = 2.43e9;        % Carrier frequency [Hz]
lambda = c / fc; % lambda [m]
halfStep = lambda / 2;


% Isolate square waveform and data
radarData = backscatter;
syncData = sync >= 0.25; % Threshold for sync signal
numSamplesPerChirp = round(chirpDuration * Fs);
numN_rp = ceil(T_rp * Fs);

% Determine the starting positions for measurements
positionChanges = diff([false; ~syncData; false]);
positionLengths = find(positionChanges < 0) - find(positionChanges > 0);
positionIndices = find(positionChanges > 0);

validPositions = positionLengths >= 2 * numN_rp; 
positionIndices = positionIndices(validPositions);
positionLengths = positionLengths(validPositions);
positionIndices = positionIndices + positionLengths - 1;
clear positionChanges validPositions positionLengths;

% Create a data matrix for positions
positionIndicesMatrix = cell2mat(arrayfun(@(x) x:x + numN_rp - 1, positionIndices(1:end-1), 'UniformOutput', false));
radarData = radarData(positionIndicesMatrix);
syncData = syncData(positionIndicesMatrix);
clear positionIndicesMatrix;

% Number of valid positions
numPositions = size(radarData, 1);

% Process each position's data
integratedData = NaN(numPositions, numSamplesPerChirp);
for posIndex = 1:numPositions
    % Identify beginnings of upchirps
    upchirpIndices = diff([false, syncData(posIndex, :), false]);
    upchirpLengths = find(upchirpIndices < 0) - find(upchirpIndices > 0);
    upchirpStartIndices = find(upchirpIndices > 0);
    upchirpStartIndices = upchirpStartIndices(upchirpLengths >= 0.9 * numSamplesPerChirp);

    % Integrate upchirps
    upchirpData = cell2mat(arrayfun(@(x) radarData(posIndex, x:x + numSamplesPerChirp - 1), upchirpStartIndices, 'UniformOutput', false));
    integratedData(posIndex, :) = mean(reshape(upchirpData', numSamplesPerChirp, length(upchirpStartIndices))', 1); %#ok<UDIM> 
end
clear upchirpStartIndices upchirpData radarData;

%%  HILBERT TRANSFORM
% Perform FFT on each row
fftData = fft(integratedData, samplesPerChirp, 2);

% Inverse FFT on positive frequencies
fftData = fftData(:, 1:end/2);
ifftData = ifft(fftData, samplesPerChirp, 2);

% Replace NaN values with a small number
ifftData(isnan(ifftData)) = 1e-30;

%%  APPLYING HANN WINDOW
numOfPositions = size(integratedData, 1);
windowedData = zeros(numOfPositions, samplesPerChirp);
hannWindow = hann(samplesPerChirp)';
for k = 1:numOfPositions
    windowedData(k, :) = ifftData(k, :) .* hannWindow;
end

% Compute wavenumber arrays
krValues = linspace(4 * pi / c * f_start, 4 * pi / c * f_stop, samplesPerChirp);
crossRangeValues = linspace(-(halfStep * numOfPositions) / 2, (halfStep * numOfPositions) / 2, numOfPositions);
figure(), clf()
imagesc(krValues, crossRangeValues, angle(windowedData))
colorbar;

%%  RANGE MIGRATION ALGORITHM
% Zero-padding
paddingSize = 2048;
paddedZeros = zeros(paddingSize, size(windowedData, 2));
for i = 1:size(windowedData, 2)
    paddingIndex = round((paddingSize - size(windowedData, 1)) / 2);
    paddedZeros(paddingIndex + 1:(paddingIndex + size(windowedData, 1)), i) = windowedData(:, i);
end
paddedWindowedData = paddedZeros;

%%  Compute the FFT of the padded data
dataMatrix = fftshift(fft(paddedWindowedData, [], 1), 1);
kxValues = linspace((-2 * pi / lambda), (2 * pi / lambda), (size(dataMatrix, 1)));

% Plot the magnitude of the cross-range FFT in dB
figure();
imagesc(krValues, kxValues, 20 * log10(abs(dataMatrix)), [max(max(20 * log10(abs(dataMatrix)))) - 40, max(max(20 * log10(abs(dataMatrix))))]);
xlabel('K_r (rad/m)'); ylabel('K_x (rad/m)');
colorbar;

% Plot the phase of the cross-range FFT
figure();
imagesc(krValues, kxValues, angle(dataMatrix));
xlabel('K_r (rad/m)'); ylabel('K_x (rad/m)');
colorbar;

% Stolt interpolation
kyValues = zeros(paddingSize, samplesPerChirp);
interpDataMatrix = zeros(paddingSize, paddingSize / 2);
for i = 1:paddingSize
    kyValues(i, :) = sqrt(krValues.^2 - kxValues(i)^2);
end

minKy = floor(min(min(kyValues)));
maxKy = ceil(max(max(kyValues)));
k_y_values = linspace(minKy, maxKy, paddingSize / 2);

for i = 1:paddingSize
    interpDataMatrix(i, :) = (interp1(kyValues(i, :), dataMatrix(i, :), k_y_values));
end

interpDataMatrix(isnan(interpDataMatrix)) = 1e-30;

% Plot the phase after Stolt interpolation
figure();
imagesc(k_y_values, kxValues, angle(interpDataMatrix)); 
xlabel('K_y (rad/m)'); ylabel('K_x (rad/m)');
colorbar;

% Perform 2D inverse FFT
ifftDataMatrix = ifft2(interpDataMatrix, (paddingSize * 4), ((paddingSize / 2) * 4));

% Calculate range and cross-range parameters
deltaFrequency = c * (maxKy - minKy) / (4 * pi);
R_max = (4 * c * paddingSize / 2) / (2 * deltaFrequency);
rail_R_max = paddingSize * halfStep;

constantFactor = 3;

d_range_1 = 1;
d_range_2 = 100;
cRangeStart = -25; 
cRangeEnd = 25;  
flippedImage = fliplr(rot90(ifftDataMatrix)); 
dIdxStart = round((size(flippedImage, 1) / R_max) * d_range_1 * constantFactor);
dIdxEnd = round((size(flippedImage, 1) / R_max) * d_range_2 * constantFactor);
cIdxStart = round((size(flippedImage, 2) / rail_R_max) * (cRangeStart + (rail_R_max / 2)));
cIdxEnd = round((size(flippedImage, 2) / rail_R_max) * (cRangeEnd + (rail_R_max / 2)));

truncatedImage = flippedImage(dIdxStart:dIdxEnd, cIdxStart:cIdxEnd);
downRangeValues = linspace(-d_range_1, -d_range_2, size(truncatedImage, 1));
crossRangeValuesFinal = linspace(cRangeStart, cRangeEnd, size(truncatedImage, 2));

% Final processing of the truncated image
finalTruncatedImage = truncatedImage;
for j = 1:size(truncatedImage, 2)
    finalTruncatedImage(:, j) = (truncatedImage(:, j)') .* (abs(downRangeValues)).^2;
end
finalTruncatedDb = 20 * log10(abs(finalTruncatedImage));

% Display the final result
figure();
imagesc(crossRangeValuesFinal, downRangeValues, finalTruncatedDb);  
xlim([-70 70]); clim([-60 -15]); colorbar;
