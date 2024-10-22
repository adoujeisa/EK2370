clear all;
close all;
clc;

% 参数设置
c = 3e8;                  % 光速 [m/s]
f_center = 2.45e9;       % 载波频率 [Hz]

% 假设相位数据是从你的 USRP 数据中计算得到的
% phase_data = angle(demodulated_data);  % 示例：解调数据的相位

% 这里以随机生成的相位数据作为示例
Fs = 1000;  % 采样频率
t = 0:1/Fs:20;  % 时间向量 [0, 20] 秒
phase_data = sin(2 * pi * 0.1 * t);  % 示例：假设相位数据是正弦波

% 计算范围
range = (c * phase_data) / (4 * pi * f_center);

% 绘制范围变化图
figure();
plot(t, range);
xlabel('Time [s]');
ylabel('Range [m]');
title('Respiratory Rate');
xlim([0, 20]); % 设置 x 轴范围
grid on;


