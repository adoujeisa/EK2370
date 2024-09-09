clear all;
close all;
clc;

% read the file
[cw_data,fs] = audioread('Velocity_Test_File.m4a');
% Parameters  数据反向？？？
N = length(cw_data);  % 样本数


% parameter
c = 3e8;  % 光速
Tp_values = [1e-3, 5e-3, 10e-3]; % Different pulse durations

% 遍历不同的 Tp 值，计算并绘制不同的速度谱
for i = 1:length(Tp_values)
    Tp = Tp_values(i);  % 当前的脉冲持续时间
    N = round(Tp * fs);  % 对应的样本数
    
    % 确保数据长度不超过信号总长度
    if N > length(cw_data)
        N = length(cw_data);
    end
    % 取前N个数据进行处理
    data_segment = cw_data(1:N);  
    
    % 执行FFT来获得速度谱
    doppler_spectrum = fftshift(fft(data_segment, N));
    freq_axis = linspace(-fs/2, fs/2, N);  % 频率轴
    
    % 绘制结果
    subplot(length(Tp_values), 1, i);
    plot(freq_axis, abs(doppler_spectrum));
    title(['CW Velocity Spectrum with Tp = ', num2str(Tp), ' s']);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
end


% MS clutter rejection
% 数据分帧处理
num_frames = floor(length(cw_data) / N);
frames = reshape(cw_data(1:num_frames * N), [N, num_frames]);

% 计算每帧的FFT
spectra = abs(fftshift(fft(frames, N, 1)));

% 计算静态背景（简单平均）
background = mean(spectra, 2);

% 计算移动目标信号（每帧减去静态背景）
target_spectra = spectra - background;

% 绘制结果
figure;
subplot(3,1,1);
plot(linspace(-fs/2, fs/2, N), mean(spectra, 2));
title('Average Spectrum Without MS Filtering');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

subplot(3,1,2);
plot(linspace(-fs/2, fs/2, N), mean(target_spectra, 2));
title('Average Spectrum With MS Filtering');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% 绘制一个样本帧的频谱
subplot(3,1,3);
plot(linspace(-fs/2, fs/2, N), spectra(:,1));
title('Spectrum of the First Frame');
xlabel('Frequency (Hz)');
ylabel('Amplitude');


figure;
%零填充 
% 设置零填充点数
N_fft = 2048000;  % 比信号长度大的 FFT 点数，MATLAB 会自动进行零填充

for i = 1:length(Tp_values)
    Tp = Tp_values(i);  % 当前的脉冲持续时间
    N = round(Tp * fs);  % 对应的样本数
    
    % 确保数据长度不超过信号总长度
    if N > length(cw_data)
        N = length(cw_data);
    end
    
    % 取前N个数据进行处理
    data_segment = cw_data(1:N);  
    
    % 执行FFT，MATLAB 会自动零填充到 N_fft 点
    doppler_spectrum = fftshift(fft(data_segment, N_fft));  % 使用 N_fft 执行 FFT
    freq_axis = linspace(-fs/2, fs/2, N_fft);  % 频率轴
    
   
    % 绘制结果
    subplot(length(Tp_values), 1, i);
    plot(freq_axis, abs(doppler_spectrum));
    title(['CW Velocity Spectrum ', num2str(Tp), ' s with MATLAB Auto Zero Padding']);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
end




% normalization 1

% normalization 2

% plot