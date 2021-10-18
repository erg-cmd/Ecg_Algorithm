%% Grafica de TF en Z
%
%
%


fs=240
ts = 1/fs
syms z
num=[1 0  1];
den=[3 0 -1];
Hz=tf2zp(num,den);
figure(1)
zplane(num,den,'r')
title('Polos y ceros de Hz=(z^2+1) /(3 z^2 - 1)')
fvtool(num,den)