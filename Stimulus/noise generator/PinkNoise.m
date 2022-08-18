clear all
close all
Data=zeros(10,10)
for i=1:10
 filei=['Noise' num2str(i) '.wav'];   
 fs = 44.1e3;
 duration = 0.6;
 Octa=["1 octave" "2/3 octave" "1/2 octave" "1/3 octave" "1/6 octave" "1/12 octave" "1/24 octave" "1/48 octave"];
 Roct=randi(length(Octa));
 oscillationFrequency = randi([100,2000],1,1);% about amplitude of dB
 centerFreq = randi([100,4000],1,1);% center frequency of filter
 bw = Octa(Roct);% bandwidth of filter
 octFilt = octaveFilter(centerFreq,bw,'SampleRate',fs);
 Data(i,1)=i;
 Data(i,2)=oscillationFrequency;
 Data(i,3)=centerFreq;
 Data(i,4)=Roct;
 y = pinknoise(duration*fs);
 y = octFilt(y);
 y=y.*reshape(sin((1:numel(y))*2*pi/(oscillationFrequency*fs/1000)),size(y,1),size(y,2));

 audiowrite(filei,y,fs);
 %----------------------------average power spectral density-------------------
 [~,freqVec,~,psd] = spectrogram(y,round(0.05*fs),[],[],fs);
 meanPSD = mean(psd,2);
 semilogx(freqVec,db(meanPSD,"power"))
 xlabel('Frequency (Hz)')
 ylabel('PSD (dB/Hz)')
 title('Power Spectral Density of Pink Noise (Averaged)')
 grid on
end