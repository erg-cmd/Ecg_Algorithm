%%
%
%
%

y= 5 * cos(w1.*t) + 5*0.8*cos(w1.*t).*cos(w2.*t);
figure(1);
plot(t,y);
a = sqrt(68 + y.^2);
b = sqrt(12^2 + y.^2);
hold on; grid on;
plot(t,a,'r');
plot(t,b,'g');



y= 5 * cos(w1.*t) + 5*0.8*cos(w1.*t).*cos(w2.*t) + 15);
figure(2);
plot(t,y);
a = sqrt(68 + y.^2);
b = sqrt(12^2 + y.^2);
hold on; grid on;
plot(t,a,'r');
plot(t,b,'g');


y= 5 * cos(w1.*t) + 5*0.8*cos(w1.*t).*cos(w2.*t) + sin(2*pi*5.*t);
figure(3);
plot(t,y);
a = sqrt(68 + y.^2);
b = sqrt(12^2 + y.^2);
hold on; grid on;
plot(t,a,'r');
plot(t,b,'g');

y= 5 * cos(w1.*t) + 5*0.8*cos(w1.*t).*cos(w2.*t) + 3*sin(2*pi*10.*t);
figure(4);
plot(t,y);
a = sqrt(68 + y.^2);
b = sqrt(12^2 + y.^2);
hold on; grid on;
plot(t,a,'r');
plot(t,b,'g');
