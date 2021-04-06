% ����ʱ������
fs = 48000;      % [Hz] �źŲ���Ƶ��
T = 2.5;          % [s] �źų��� 
x = 0:1/fs:T;   % [s] ʱ������
% �����ź�����
f = 6000;        % [Hz] �ź�Ƶ��
l = 10000; % ���
y = 1*sin(2*pi*f*x);
% �����Ƶ�ļ�
fname = [num2str(f),'.wav'];        % �趨�ļ����� ע���ʽ
audiowrite(fname,y,fs);     % ����ļ�

[X,FS]=audioread([num2str(f), '.WAV']); % �� WAV �ļ�ת���ɱ���
t=(0:length(X)-1)/FS; % ��������ʱ��
xf=fft(X); % ������Ҷ�任��ԭƵ��
fm=10000*length(xf)/FS; % ȷ����Ƶ��ͼ������Ƶ��
f=(0:fm)*FS/length(xf); % ȷ����Ƶ��ͼ��Ƶ�ʿ̶�
% plot(f,abs(xf(1:length(f)))); % ����ԭ����Ƶ��ͼ


% y = 0.2*sin(2*pi*25*(1:12000)/200);
% fs = 1000;
% sound(y)
% fname = 'high.wav';        % �趨�ļ����� ע���ʽ
% audiowrite(fname,y,fs); 
% [X,FS]=audioread('high.WAV'); % �� WAV �ļ�ת���ɱ���
% t=(0:length(X)-1)/FS; % ��������ʱ��
% % plot(t,X); % ����ԭ����ͼ
% xf=fft(X); % ������Ҷ�任��ԭƵ��
% plot(f,abs(xf(1:length(f)))); % ����ԭ����Ƶ��ͼ



% %���� FDATool ���һ�� LowpassButterworth �˲���
% %ָ�� FS=22050Hz Fp=1000Hz Ap=1dB Fs=3000Hz As=20dB
% B =[0.0062,0.0187,0.0187,0.0062]; % ����ϵ��
% A =[1,-2.1706,1.6517,-0.4312]; % ��ĸϵ��
% Y=filter(B,A,X); % ʵ�������˲�
% t=(0:length(X)-1)/FS; % ��������ʱ��
% subplot(2,2,1);plot(t,X); % ����ԭ����ͼ
% title(' ԭ�źŲ���ͼ '); % �ӱ���
% subplot(2,2,3);plot(t,Y); % �����˲�����ͼ
% title(' �˲�����ͼ '); % �ӱ���
% xf=fft(X); % ������Ҷ�任��ԭƵ��
% yf=fft(Y); % ������Ҷ�任���˲���Ƶ��
% fm=3000*length(xf)/FS; % ȷ����Ƶ��ͼ������Ƶ��
% f=(0:fm)*FS/length(xf); % ȷ����Ƶ��ͼ��Ƶ�ʿ̶�
% subplot(2,2,2);
% plot(f,abs(xf(1:length(f)))); % ����ԭ����Ƶ��ͼ
% title(' ԭ�ź�Ƶ��ͼ '); % �ӱ���
% subplot(2,2,4);plot(f,abs(yf(1:length(f)))); % �����˲���Ƶ��ͼ
% title(' �˲����ź�Ƶ��ͼ '); % �ӱ���