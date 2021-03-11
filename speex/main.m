clc;clear ;

[rrin,fs] = audioread('ref16k.wav');
[ssin,fs] = audioread('mic16k.wav');

Fs =fs;
u = rrin ;
d = ssin ;
filter_length = 4096 ;
frame_size = 128;

speex_mdf_out = speex_mdf(Fs, u, d, filter_length, frame_size, 'dbg_var_name');
 
 
 