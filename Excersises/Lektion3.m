%% Lektion 3

%% 2.25

N = 60;
p = 2; %length after decimation
n = 0:N-1;
w0 = 2*pi/5;
w = (-N/2:N/2-1)*2*pi/N;

x = sin(w0*n);
X = fftshift(fft(x));

x2 = x(1:p:N);
X2 = fftshift(fft(x2));

M=length(x2);
w2 =( -M/2:M/2-1)*2*pi/M;

figure(2)
subplot(2,1,1);
stem(w,abs(X))
subplot(2,1,2);
stem(w2,abs(X2))


%% 3.2
load sig30

y1 = [y' zeros(1,10000)];
Y1 = fftshift(fft(y1));
N =length(Y1);

w = (-N/2:N/2-1)*2*pi/N

figure(3)
plot(w,abs(Y1))



%% 4.16

load sig40

T = 1;
N = length(s);

w0 = 0.3;
wn = w0*T/pi;

[B A] = butter(10,wn,'low')

s1a = filter(B,A,s);
s1b = filtfilt(B,A,s);
n = (0:N-1)

figure(4)
subplot(3,1,1)
plot(n,s1)
subplot(3,1,2)
plot(n,s1a)
subplot(3,1,3)
plot(n,s1b)

[B A] = butter(10,wn,'high')

s2a = filter(B,A,s);
s2b = filtfilt(B,A,s);
n = (0:N-1)

figure(5)
subplot(3,1,1)
plot(n,s2)
subplot(3,1,2)
plot(n,s2a)
subplot(3,1,3)
plot(n,s2b)

%% 4.18