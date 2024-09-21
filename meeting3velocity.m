clear all;
close all;
clc;

%% COMPUTE VELOCITY FROM AUDIO FILE
% Read the audiofile 
[y, fs] = audioread('Recording (4).m4a'); 

% Take care of data inversion by the sound card
data = -y(:, 1); % Intensity of the received signal

% Parameters
c = 3e8;                % Speed of light [m/s]
f_center = 2.43e9;            % Center Frequency [Hz]
Tp = 0.1;                     % Pulse width [s]
N = round(Tp * fs);           % Number of samples per pulse
N_total = length(data);       % Total number of samples in read file
T = N_total / fs;             % Total time duration of signals in read file
num_pulse = round(T / Tp);    % Total number of pulses

% Initialize time-domain matrix
mat_time = zeros(num_pulse, N); 
for i = 1:num_pulse-1
    mat_time(i, :) = data(1 + N * (i-1) : i * N); % Assign amplitude
end

%% Mean subtraction clutter removal (MS)
data_mean = mean(data); % Mean of the entire data
mat_time = mat_time - data_mean; % Remove clutter by subtracting the mean
%%


%%zero padding for fft 使用 reshape 函数将原始数据矩阵 data 重塑为时间域矩阵 mat_time，每一行对应一个脉冲。
for i = 1:num_pulse-1
    mat_freq(i,:)= 20*log(abs((fft(mat_time(i,:),6*N))));
 end


%% FFT and Velocity Calculation
fft_data = fft(mat_time, [], 2); % Apply FFT along the second dimension (pulses)
fft_data_db = 20 * log10(abs(fft_data)); % Convert magnitude to dB

 %% Normalize 1
fft_data_norm = fft_data_db - max(max(fft_data_db));
 %% Normalize 2
% Calculate the maximum value for each column
max_vals = max(fft_data_db, [], 2);
% Subtract each column's maximum value from each element in that column
fft_data_norm_2 = fft_data_db - max_vals;


% Velocity array
f_doppler = linspace(0, fs/2, 6*N); % Doppler frequency axis
velocity = f_doppler * c / (2 * f_center); % Convert Doppler frequency to velocity

% Time array
time_array(:,1) = Tp:Tp:Tp*num_pulse;
%time_array = (0:num_pulse-2) * Tp; % Time array for the y-axis

%% Plot using imagesc
figure (1);
imagesc(velocity, time_array, fft_data_norm); % Plot the normalized FFT data
xlabel('Velocity (m/s)');
ylabel('Time (s)');
title('Velocity-Time Map with pulse time = 0.1s noormalization 1 with 6N');
colorbar; % Add color bar to represent intensity

xlim([0 20]); % Set the range for the x-axis (velocity) from 0 to 30 m/s
clim([-30 0]); % Set the color axis range for the colorbar, e.g., from -45 dB to 0 dB

figure (2);
imagesc(velocity, time_array, fft_data_norm_2); % Plot the normalized FFT data
xlabel('Velocity (m/s)');
ylabel('Time (s)');
title('Velocity-Time Map with pulse time = 0.1s noormalization 2 with 6N');
colorbar; % Add color bar to represent intensity

xlim([0 20]); % Set the range for the x-axis (velocity) from 0 to 30 m/s
clim([-10 0]); % Set the color axis range for the colorbar, e.g., from -15 dB to 0 dB


% 捕获数据图
figure(3);
t = (0:N_total-1) / fs;  % 创建时间轴
plot(t, data);
xlabel('Time (s)');
ylabel('Magnitude of the sigal');
title('captured data plot');
grid on;

% 取最强散射体
[max_val, max_idx] = max(fft_data_db, [], 2); % 找到每一帧中最大的频谱分量
strongest_velocity = velocity(max_idx); % 取出对应的速度

% 速度-时间图
figure(4);
plot(time_array, strongest_velocity);
xlabel('Time (s)');
ylabel('V (m/s)');
title('v-t');
grid on;

% CW Movie of real-time/pseudo real-time velocity measurements
% 为什么这俩的scale不一样啊？？？？？？？
video_cw = VideoWriter('CW_Velocity_Measurement.avi');
open(video_cw);

figure;
velocity_trace = []; % 用于存储所有的速度点
time_trace = [];     % 用于存储所有的时间点

for i = 1:num_pulse-1
    % 计算当前脉冲的FFT
    fft_current = fft(mat_time(i,:), 6*N);
    fft_current_db = 20 * log10(abs(fft_current));
    
    % 找出最大散射体的速度
    [~, max_idx] = max(fft_current_db);
    current_velocity = velocity(max_idx);
    
    % 记录当前时间和速度
    velocity_trace = [velocity_trace, current_velocity];
    time_trace = [time_trace, Tp * i];
    
    % 绘制当前所有的点，而不是逐个点
    plot(time_trace, velocity_trace, 'ro-');
    xlabel('time (s)');
    ylabel('v (m/s)');
    title('real-time/pseudo real-time velocity measurements');
    xlim([0, num_pulse * Tp]);
    ylim([0, 2]);  % 调整Y轴范围，使其更符合实际速度
    grid on;
    
    % 将当前帧写入视频
    drawnow; % 确保图形刷新
    frame = getframe(gcf); % 获取完整的图像
    writeVideo(video_cw, frame);
end
close(video_cw);
