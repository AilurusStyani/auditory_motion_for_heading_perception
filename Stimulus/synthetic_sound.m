% 生成时间序列
fs = 48000;      % [Hz] 信号采样频率
T = 2.5;          % [s] 信号长度 
x = 0:1/fs:T;   % [s] 时间序列
% 生成信号序列
f = 6000;        % [Hz] 信号频率
l = 10000; % 振幅
y = 1*sin(2*pi*f*x);
% 输出音频文件
fname = [num2str(f),'.wav'];        % 设定文件名称 注意格式
audiowrite(fname,y,fs);     % 输出文件

[X,FS]=audioread([num2str(f), '.WAV']); % 将 WAV 文件转换成变量
t=(0:length(X)-1)/FS; % 计算数据时刻
xf=fft(X); % 作傅里叶变换求原频谱
fm=10000*length(xf)/FS; % 确定绘频谱图的上限频率
f=(0:fm)*FS/length(xf); % 确定绘频谱图的频率刻度
% plot(f,abs(xf(1:length(f)))); % 绘制原波形频谱图


% y = 0.2*sin(2*pi*25*(1:12000)/200);
% fs = 1000;
% sound(y)
% fname = 'high.wav';        % 设定文件名称 注意格式
% audiowrite(fname,y,fs); 
% [X,FS]=audioread('high.WAV'); % 将 WAV 文件转换成变量
% t=(0:length(X)-1)/FS; % 计算数据时刻
% % plot(t,X); % 绘制原波形图
% xf=fft(X); % 作傅里叶变换求原频谱
% plot(f,abs(xf(1:length(f)))); % 绘制原波形频谱图



% %利用 FDATool 设计一个 LowpassButterworth 滤波器
% %指标 FS=22050Hz Fp=1000Hz Ap=1dB Fs=3000Hz As=20dB
% B =[0.0062,0.0187,0.0187,0.0062]; % 分子系数
% A =[1,-2.1706,1.6517,-0.4312]; % 分母系数
% Y=filter(B,A,X); % 实现数字滤波
% t=(0:length(X)-1)/FS; % 计算数据时刻
% subplot(2,2,1);plot(t,X); % 绘制原波形图
% title(' 原信号波形图 '); % 加标题
% subplot(2,2,3);plot(t,Y); % 绘制滤波波形图
% title(' 滤波后波形图 '); % 加标题
% xf=fft(X); % 作傅里叶变换求原频谱
% yf=fft(Y); % 作傅里叶变换求滤波后频谱
% fm=3000*length(xf)/FS; % 确定绘频谱图的上限频率
% f=(0:fm)*FS/length(xf); % 确定绘频谱图的频率刻度
% subplot(2,2,2);
% plot(f,abs(xf(1:length(f)))); % 绘制原波形频谱图
% title(' 原信号频谱图 '); % 加标题
% subplot(2,2,4);plot(f,abs(yf(1:length(f)))); % 绘制滤波后频谱图
% title(' 滤波后信号频谱图 '); % 加标题