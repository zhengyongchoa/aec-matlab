% 
% author: zhengyongchao 2018/2/27
%   Jean-Marc Valin, Member, IEEE
%  On Adjusting the Learning Rate in Frequency Domain
%   Echo Cancellation With Double-Talk  


clear ; clc; close all;
 %% 
fle_mic = audioread('mic16k.wav');
fle_ref =  audioread('ref16k.wav');

  %%
  Fs = 16000;
  x = fle_ref;
  y = fle_mic;  

 %%
N=512; M=64; P=N/M;
yeo = y;
e=zeros(size(y));
 %%
wF=zeros(2*M,P);
xF=zeros(2*M,P);
mu=0.1;
for n=M+1:M:length(x)-M
    xF=[fft(x(n-M:n+M-1)) xF(:,[1:end-1])];    %第一列更新fftx，其他列后移
    yhat=ifft(sum((wF.*xF).').');              %计算滤波器与信号的频域乘法（按行求和）
    yhat=real(yhat(M+1:end));                  %舍弃圆周卷积的数据
    E=yeo(n:n+M-1)-yhat;                       %误差信号
    e(n:n+M-1)=E;
    MU=mu*(sum((abs(xF).^2)')'+0.1).^(-1);        %归一化步长
    EF=fft([zeros(M,1); E]);                      %补零凑长度
    wF=wF+diag(MU.*EF)*conj(xF);                  %更新滤波系数
    waux=real(ifft(wF)); wF=fft([waux(1:M,:); zeros(M,P)]);%正反变换补零


end
 %%
 audiowrite('Residual sound.wav',e,Fs);
 audiowrite('mic sound.wav',y,Fs);


