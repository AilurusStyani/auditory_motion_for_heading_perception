clear all
close all
Data=zeros(100,100);
%for i=1:8
       y=zeros(100,100);
       %filei=['Noise' num2str(1000*(i-1)+100*(j-1)+l) '.wav'];   
       fs = 44.1e3;
<<<<<<< Updated upstream
       duration = 0.6;
=======
       duration = 1.2;
>>>>>>> Stashed changes
       Octa=["1 octave" "2/3 octave" "1/2 octave" "1/3 octave" "1/6 octave" "1/12 octave" "1/24 octave" "1/48 octave"];
       Roct=randperm(length(Octa));
       %oscillationFrequency = randi([100,900],1,1);% about amplitude of dB, recommend range 100 900
       oscillationFrequency = [20,50,100,200,300,400,500,600];
       osiFreGroup=oscillationFrequency(randperm(length(oscillationFrequency)));
       %centerFreq = randi([1000,10000],1,1);% center frequency of
       %filter,recommend range 1000~20000
       centerFreq = [500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000];
       cenFreGroup=centerFreq(randperm(length(centerFreq)));
       cenFreGroup=cenFreGroup(1:length(Octa));
       bw = Octa(Roct);% bandwidth of filter
      for i=1:8
       octFilt = octaveFilter(cenFreGroup(i),bw(i),'SampleRate',fs);
% %        Data={100*i+10*j+l, oscillationFrequency(j), centerFreq(l), Roct(i)}l
% %        Data(l+,1)=100*i+10*j+l;
% %        Data(i,2)=oscillationFrequency(j);
% %        Data(i,3)=centerFreq(l);
% %        Data(i,4)=Roct(i);
       y = pinknoise(duration*fs);
       y= octFilt(y);
       y=y.*reshape(sin((1:numel(y))*2*pi/(osiFreGroup(i)*fs/1000)),size(y,1),size(y,2));
       plot(y);
<<<<<<< Updated upstream
       filei=['D:\LQY\auditory_motion_for_heading_perception220215\Stimulus\220812 05\Noise' num2str(i+56 ) '.wav'];
=======
       filei=['D:\LQY\auditory_motion_for_heading_perception220215\Stimulus Lifetime\1200ms\Noise' num2str(i+80) '.wav'];
>>>>>>> Stashed changes
       audiowrite(filei,y,fs);
      end
%  ----------------------------average power spectral density-------------------
%   [~,freqVec,~,psd] = spectrogram(y,round(0.05*fs),[],[],fs);
%   meanPSD = mean(psd,2);
%   semilogx(freqVec,db(meanPSD,"power"))
%   xlabel('Frequency (Hz)')
%   ylabel('PSD (dB/Hz)')
%   title('Power Spectral Density of Pink Noise (Averaged)')
%   grid on

%end