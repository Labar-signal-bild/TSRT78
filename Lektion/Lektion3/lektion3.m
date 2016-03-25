%% TSRT78 Lektion 3
initcourse('TSRT78')
%% 2.25
T = 1;
w0 = 2*pi/5;
N = 60;
n = (0:N-1)*T;

xn = sin(w0*n);
w = (-N/2:(N/2-1))*2*pi/(N*T);
Xn = fftshift(fft(xn)*T);

% a) Decimate the signal a factor of two
p = 2;

xhat = xn(1:p:length(xn));
w1 = (-N/p/2:(N/p/2-1))*2*pi/(N/p*T);
Xhat = fftshift(fft(xhat)*p*T);


% b) Decimate a factor of 3
p = 3;

xhat2 = xn(1:p:length(xn));
w2 = (-N/p/2:(N/p/2-1))*2*pi/(N/p*T);
Xhat2 = fftshift(fft(xhat2)*p*T);


figure(1)
subplot(1,2,1)
stem(w, abs(Xn))
hold on
subplot(1,2,1)
stem(w1, abs(Xhat))

subplot(1,2,2)
stem(w, abs(Xn))
hold on
subplot(1,2,2)
stem(w2, abs(Xhat2))


% c) Decimate a factor of 7(!!)













