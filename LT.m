%% Prueba de la length Transform
%
%
%
%
A = 5;
f = 50;
w = 2*pi*f;
t = 0:0.01:2*pi;
y = A*(sind(w*t)+0.38*sind(w.*t).*sind(10.*w.*t))+5;
LTwindow = ceil(0.13 * f);
ltfs = 500000*A*A/f;
lt = zeros(length(y));			%Vector "lt" del tamaño de muestras llenado con 0

for m = 1: length(y)-LTwindow + 1
	dy = sum(y(m:m+LTwindow-1));
	lt(m+LTwindow-1) = sqrt(ltfs+power(dy,2));  % hacemos la length transform
end

figure(1);grid;
plot(t,y,'r');
hold on;
plot(t,lt,'g');


