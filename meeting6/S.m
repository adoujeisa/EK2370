clear all;
close all;
clc;


% 文件夹路径及文件名
folder_path = 'D:\EK2370\sfile'; % 修改为实际路径
file_names = {'amplifier1.s2p', 'amplifier2.s2p', 'amplifier3.s2p', 'amplifier4.s2p', 'amplifier5.s2p', 'amplifier6.s2p', 'amplifier7.s2p'};

% 初始化存储数据的数组
s11_data = {};
s21_data = {};
s12_data = {};
frequency_data = {};

% 读取每个 .s2p 文件并提取S参数
for i = 1:length(file_names)
    file_path = fullfile(folder_path, file_names{i});
    s_params = sparameters(file_path);
    
    % 提取频率和S参数数据
    frequency_data{i} = s_params.Frequencies; % 提取频率 [Hz]
    s11_data{i} = 20*log10(abs(squeeze(s_params.Parameters(1,1,:)))); % 提取并转换为dB
    s21_data{i} = 20*log10(abs(squeeze(s_params.Parameters(2,1,:))));
    s12_data{i} = 20*log10(abs(squeeze(s_params.Parameters(1,2,:))));
end

% 绘制S11, S21, S12
figure;
for i = 1:length(file_names)
    % 绘制S11
    subplot(3, 1, 1);
    plot(frequency_data{i}/1e9, s11_data{i}, 'DisplayName', ['File ' num2str(i)]);
    hold on;
    title('S11 vs Frequency');
    xlabel('Frequency (GHz)');
    ylabel('S11 (dB)');
    grid on;
    legend('show');
    
    % 绘制S21
    subplot(3, 1, 2);
    plot(frequency_data{i}/1e9, s21_data{i}, 'DisplayName', ['File ' num2str(i)]);
    hold on;
    title('S21 vs Frequency');
    xlabel('Frequency (GHz)');
    ylabel('S21 (dB)');
    grid on;
    legend('show');
    
    % 绘制S12
    subplot(3, 1, 3);
    plot(frequency_data{i}/1e9, s12_data{i}, 'DisplayName', ['File ' num2str(i)]);
    hold on;
    title('S12 vs Frequency');
    xlabel('Frequency (GHz)');
    ylabel('S12 (dB)');
    grid on;
    legend('show');
end
