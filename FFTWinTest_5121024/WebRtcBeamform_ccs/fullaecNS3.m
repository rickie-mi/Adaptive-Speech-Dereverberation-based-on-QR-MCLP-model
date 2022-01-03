close all
clear all
clc
% [data1,fs1]=audioread('ch1_out.wav');
% data1=data1*32768;
fid=fopen('sweep20_20K_30s_out.pcm', 'rb');      %替换dsp的fft后的文件
data1=fread(fid,inf,'int16');
data1=data1(513:end);
fclose(fid); 
% [data2,fs2]=audioread('ch1_out.wav');
% % data1=data1;
fid=fopen('sweep20_20K_30s.pcm', 'rb');
data2=fread(fid,inf,'int16');
fclose(fid); 
len=min(length(data1),length(data2));
plot(data1(1:len)-data2(1:len))

% figure
% subplot(211)
% plot(data1,'r')
% hold on
% plot(data2,'k')




