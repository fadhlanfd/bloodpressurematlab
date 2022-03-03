clc
clear
close all
video_file= 'rekamantekanandarah_4G9D54ZT.mp4';
[y,fps] = acquire(video_file)
save('data_test.mat','y','fps');
process;