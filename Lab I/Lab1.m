addpath Sound;

%% 4.2 Wistle

[whistle, fSamp] = audioread('whistle4.wav'); % extract data and sampling frequency
%%
x = whistle(8000*5+1:8000*7); % extract the two best seconds BETTER THEN 3-5????

%%
nSamp = size(x,1); % number of samples
t = (0:nSamp-1)/fSamp; % time vector in seconds

figure(1);clf(); 
plot(t, x)
xlabel('time in seconds')
ylabel('recorded signal')
%%
N = 16000;
X = fft(x);
Xsamp = length(X);

% 1: Energy time domain 
Etot1 = sum(abs(x).^2); % Definition of the energy in the time domain

[B,A] = butter(5,[(1892-30)/8000 (1892+30)/8000]);

x_t = filter(B,A,x);

Edom1 = sum(abs(x_t).^2);

figure(10)
plot(1:N,abs(fft(x_t)))
figure(11)
plot(1:N,abs(fft(x)))



% 2: Energy frequency domain
X_lowMax = find(abs(X(1:Xsamp/2)) == max(abs(X(1:Xsamp/2))));
X_highMax = find(abs(X(Xsamp/2+1:Xsamp)) == max(abs(X(Xsamp/2+1:Xsamp))))+Xsamp/2;

band_pass = zeros(1, N); % Construct a band-pass filter to obtain the dominant frequences
band_pass(X_lowMax-30:X_lowMax+30) = 1;
band_pass(X_highMax-30:X_highMax+30) = 1;

domsamples = band_pass'.*X;


Etot2 = 1/N*sum(abs(X).^2);
Edom2 = 1/N*sum(abs(domsamples).^2);

%%
sound(abs(domsamplestime),8000); 
%% 3: Harmonic distortion 

harmonic1 = 1-Edom1/Etot1;
harmonic2 = 1-Edom2/Etot2;

%% 4: AR


mo = ar(x,2); % AR model of order 2
A = mo.a;
[R,P,K] = residue(1,A); 

placement = norm(P(1)); % The length from the origin to the pole
distance = 1-placement; % The distance from the unit circle to the pole


%% 5: Estimate the dominant frequancy paratertic/non-parametric

% One of them should be a bode plot. Which one?? Both?

spect_val = etfe(x, 100, [], 1/8000); % The non-parametric method

figure(1); bode(mo, 'r', spect_val, 'b');
legend('AR-model', 'Original Signal');



%%
%------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------

%% Vowel
%-------------------------------------------------------------------------------
%----------------------AAA------------------------------------------------------
%-------------------------------------------------------------------------------



[aaa, fSamp] = audioread('A2.wav'); % extract data and sampling frequency
%sound(aaa,fSamp); % make sure that the quality of the recording is okay
nSamp = size(aaa,1); % number of samples
t = (0:nSamp-1)/fSamp; % time vector in seconds

figure(1);clf();
plot(t, aaa)
xlabel('time in seconds')
ylabel('recorded signal') % axis description is important!

aa = aaa(floor(8000*3.27)+1:floor(8000*5.27));

nSamp = size(aa,1); % number of samples
t = (0:nSamp-1)/fSamp; % time vector in seconds

figure(1);clf();
plot(t, aa)
xlabel('time in seconds')
ylabel('recorded signal') % axis description is important!

%% Cross validation

T= 1/8000;
N=length(aa);

edata = iddata(aa(1:2*floor(N/3)),[],T); % Estimated data
vdata = iddata(aa(2*floor(N/3)+1:N),[],T); % Validation data, compare w this
WNe=[];
WNv=[];
for n=1:20
mo = ar(edata,n);
WNe = [WNe mo.NoiseVariance];
res=resid(mo,vdata);
WNv = [WNv sum(res.y.*res.y)/res.N];
end
plot(1:20,WNe,1:20,WNv)

%% Spectral validation

aa_spect = etfe(vdata, 200); % Validation spectrum

mo1 = ar(edata, 1);
mo4 = ar(edata, 4);
mo7 = ar(edata, 7);
mo8 = ar(edata, 8);
mo9 = ar(edata, 9);
mo30 = ar(edata, 30);
mo20 = ar(edata, 20);

figure(2); bode(mo4, 'r', mo7, 'g', mo8, 'b', mo20, 'm', aa_spect, 'y');
legend(' = mo4', ' = mo7', ' = mo8', ' = mo20', ' = validation');

%%


mo7 = ar(edata,7);

Ts=40;
fs = 7970;
N = fs*2;
pulselength = 100;

pulsetrain = zeros(1,N);

for n=1:Ts:N
pulsetrain(n)=0.01;
end

aa_out = filter(1,mo7.a,pulsetrain);

sound(10*aa_out, 8000);

figure(1)
plot(1:length(aa),aa);%axis([5500 6000 -0.04 0.03])
hold on;
plot(1:length(aa_out),aa_out);%axis([5500 6000 -0.04 0.03])
hold off

%% O

[iii, fSamp] = audioread('I2.wav');

%sound(iii,8000)
nSamp = size(iii,1); % number of samples
t = (0:nSamp-1)/fSamp; % time vector in seconds

ii = iii(floor(8000*2.28)+1:floor(8000*4.28));

nSamp = size(ii,1); % number of samples
t = (0:nSamp-1)/fSamp; % time vector in seconds

figure(3);clf();
plot(t, ii)
xlabel('time in seconds')
ylabel('recorded signal') % axis description is important!


%% Cross validation

T= 1/8000;
N=length(ii);

edata = iddata(ii(1:2*floor(N/3)),[],T); % Estimated data
vdata = iddata(ii(2*floor(N/3)+1:N),[],T); % Validation data, compare w this
WNe=[];
WNv=[];
for n=1:20
mo = ar(edata,n);
WNe = [WNe mo.NoiseVariance];
res=resid(mo,vdata);
WNv = [WNv sum(res.y.*res.y)/res.N];
end
plot(1:20,WNe,1:20,WNv)

%%

mo20 = ar(edata,20);

ii_out = filter(1,mo20.a,pulsetrain);

sound(50*ii_out,8000);

figure(1)
plot(1:length(ii),ii);%axis([5500 6000 -0.04 0.03])
hold on;
plot(1:length(ii_out),ii_out);%axis([5500 6000 -0.04 0.03])
hold off


%------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------

%% Speech encoding as in GSM


[sen, fSamp] = audioread('Sentance5.wav'); % extract data and sampling frequency
%sound(aaa,fSamp); % make sure that the quality of the recording is okay
nSamp = size(sen,1); % number of samples
t = (0:nSamp-1)/fSamp; % time vector in seconds

sent = sen(floor(8000*1.10)+1:floor(8000*5.10));

nSamp = size(sent,1); % number of samples
t = (0:nSamp-1)/fSamp; % time vector in seconds
%sound(10*sent,8000);

%% 
order = 8;
sent_div=zeros(200,160);
for i=1:200
   sent_div(i,:)=sent(160*(i-1)+1:160*(i));
end

sent_ar = zeros(200, order+1);
for i=1:200
   temp = ar(sent_div(i,:), order);
   temp_a = poole_stable(temp.a);
   sent_ar(i,:) = temp_a; 
end

%%

yhat = zeros(200,160);

for i=1:200
    e_vec = filter(sent_ar(i,:),1,sent_div(i,:)); % m1 <-> AR model of the segment
    r = covf(e_vec',100);
    [A,D] = max(r(20:end));
    D = D+20;
    ehat = zeros(1,160);
    ehat(1:D:end)=sqrt(A);
    yhat(i,:) = filter(1,sent_ar(i,:),ehat);
end

y = zeros(1, 32000);
for i=1:200
    y(160*(i-1)+1:160*(i)) = detrend(yhat(i,:));
end

y=y./(mean(abs(y)));

figure(4); plot(y);
mean(abs(y))
%sound(y,8000)
%%
Y = fft(y);

%saturated = find(abs(Y) > 4000);

%Y(saturated) = Y(saturated)./10;

y_new = ifft(Y);

sound(y_new,8000);

%%
figure(5)
plot(1:length(Y),abs(Y));

Sent = fft(sent);
figure(6)
plot(1:length(Sent),abs(Sent));



%%
figure(3)
plot(1:32000,sent)


%TODO: Are the poles stable?? If not, mirror then in the unit circle by the
% help of poole_stable.m (make poole_stable.m)
%TODO: detend! Maybe in the vowel part as well??



