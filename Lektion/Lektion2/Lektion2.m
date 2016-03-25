%% TSRT78 Lektion 2
initcourse('TSRT78')

%% 2.3

T = 1;
N = 32;
w1 = 1;
n = (0:N-1)*T;
t = (-1000:1000)*T;

xn = cos(w1*t);
[X, theta1] = dtft(xn);
yn = [zeros(1, 100) cos(w1*n) zeros(1,100)];
[Y, theta] = dtft(yn);

%Plots of X and Y. In Y we can see leakage. The leakage is there because we
%truncate the sigal with a window.
figure(1)
subplot(2, 1, 1)
plot(theta1, abs(X));
subplot(2, 1, 2)
plot(theta, abs(Y));


%% 2.5

load power;


T = 39*10^-5;

%Downsample u2
y = u2(1:3:length(u2));

figure(2)
subplot(1, 2, 1)
plot(u2);
subplot(1, 2, 2)
plot(y);

[U2, theta1] = dtft(u2,T);
[Y, theta2] = dtft(y,3*T);

figure(3)
subplot(1, 2, 1)
plot(theta1, abs(U2)); title('U2');
subplot(1, 2, 2)
plot(theta2, abs(Y)); title('Y');

%Low pass filter u2 before downsampling

R = 3; % The down sample factor
y2 = decimate(u2, R); % Decimate ow pass filter the signal and then down sample it

[Y2, theta3] = dtft(y2,3*T);

figure(4)
subplot(1, 2, 1)
plot(theta1, abs(U2)); title('U2');
subplot(1, 2, 2)
plot(theta3, abs(Y2)); title('Y2');

% When filtering before downsampeling the signal looks better. There's not 
% as much signal outside of the Nyqvist frequency that will be alaised

%% 2.10

T = 1;
N = 16;
n = (0:N-1);

x0 = cos((2*pi/8)*n);
x1 = cos((2*pi/7)*n);

w = (0:N-1)*2*pi/(N*T);

X0 = fft(x0)*T; %Creating the DFT of the signals
X1 = fft(x1)*T;

[X0DF, W] = dtft(x0);
[X1DF, W] = dtft(x1);


figure(5)
subplot(2, 2, 1)
plot(w, abs(X0)); title('X0 and XODF');
hold on 
subplot(2, 2, 1)
plot(W, abs(X0DF)); 

subplot(2, 2, 2)
plot(w, abs(X1)); title('X1 and X1DF');
hold on
subplot(2, 2, 2)
plot(W, abs(X1DF)); 

% Recreate the DTFT using DFT by zeropadding the signal

p = 100;
x02 = [x0 zeros(1,N*(p-1))];
x12 = [x1 zeros(1,N*(p-1))];

X02 = fft(x02)*T;
X12 = fft(x12)*T;

w2 = (0:N*p-1)*2*pi/(N*p*T);

figure(5)
subplot(2, 2, 3)
plot(w2, abs(X02)); title('X02');
subplot(2, 2, 4)
plot(w2, abs(X12)); title('X12');

%% 2.22

T = 1;
N = 16;
n = (0:N-1);
w0 = 2*pi/sqrt(17);

xn = sin(w0*n);
theta = (0:N-1)*2*pi/(N*T);

X = fft(xn);

figure(6)
subplot(2, 2, 1)
plot(theta, abs(X)); title('DFT x[n]');

p1 = 16;
p2 = 64-16;
p3 = 256-16;

xn1 = [xn zeros(1, p1)];
xn2 = [xn zeros(1, p2)];
xn3 = [xn zeros(1, p3)];

X1 = fft(xn1);
X2 = fft(xn2);
X3 = fft(xn3);

theta1 = (0:N+p1-1)*2*pi/(N+p1*T);
theta2 = (0:N+p2-1)*2*pi/(N+p2*T);
theta3 = (0:N+p3-1)*2*pi/(N+p3*T);


subplot(2, 2, 2)
plot(theta1, abs(X1)); title('DFT x1[n]');
hold on 
plot(theta, abs(X));

subplot(2, 2, 3)
plot(theta2, abs(X2)); title('DFT x2[n]');
hold on 
plot(theta, abs(X));

subplot(2, 2, 4)
plot(theta3, abs(X3)); title('DFT x4[n]');
hold on 
plot(theta, abs(X));




