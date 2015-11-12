%% Lektion 2

initcourse TSRT78



%% 2.3
% b
N = 32; 
T = 1;
n = (0:N-1).*T;
w1 = 1;

xn      = cos(w1*n);
yn      = [zeros(1,64) xn zeros(1,64)];
[Y,theta]       = dtft(yn);


figure(1) 
plot(theta,abs(Y))
% c

%% 2.5
% a
load power
N = length(u2);
T = 3.9*10^-4
y = u2(1:3:u2end);
[Y thetay] = dtft(y,3*T);
[U2 thetaU2] = dtft(u2,T);

figure(2) 
subplot(2,1,1)
plot(thetay,abs(Y));title('Y');
subplot(2,1,2)
plot(thetaU2,abs(U2));title('U2')

% b
u2Lp = decimate(u2,3);

[U2Lp thetaLP] = dtft(u2Lp,3*T);

figure(3) 
subplot(2,1,1)
plot(thetaLP,abs(U2Lp));title('U2Lp');
subplot(2,1,2)
plot(thetaU2,abs(U2));title('U2');

% much smother then before the Lp and we don't get "false" frequencies in
% the middle.

%% 2.10
% a
N = 16;
n = 0:N-1;
w = (2*pi/N)*n;

x0 = cos(2*pi/8*n);
x1 = cos(2*pi/7*n);

X0F = fft(x0);
X1F = fft(x1);

w = (0:N-1)*2*pi/N;

[X0 t0]  = dtft(x0)
[X1 t1]  = dtft(x1)


figure(4)
subplot(2,1,1)
plot(w,abs(X0F)); title('X0F');
subplot(2,1,2)
plot(w,abs(X1F)); title('X1F');
subplot(2,1,3)
plot(t0,abs(X0)); title('X0');
subplot(2,1,4)
plot(t1,abs(X1)); title('X1');

% b 
p = 4;
K = N(p-1);

%x02 =
%x01 =



%% 2.22
%b

N = 16;
n =0:N-1;
w0 = 2*pi*sqrt(17);
x = sin(w0*n);

p1 = 2;
p2 = 4;
p3 = 16;
x1 = [x zeros(1,N*(p1-1))];
x2 = [x zeros(1,N*(p2-1))];
x3 = [x zeros(1,N*(p3-1))];












