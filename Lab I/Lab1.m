%% Wistle
addpath Sound;
%%
[whistle, fSamp] = audioread('whistle4.wav'); % extract data and sampling frequency
%%
x = whistle(8000*3+1:8000*5);
%%
nSamp = size(x,1); % number of samples
t = (0:nSamp-1)/fSamp; % time vector in seconds

figure(1);clf();
plot(t, x)
xlabel('time in seconds')
ylabel('recorded signal') % axis description is important!
%%
N = 16000;
X = fft(x);
Xsamp = length(X);

% E time domain
Etot1 = sum(abs(x).^2);

X_lowMax = find(abs(X(1:Xsamp/2)) == max(abs(X(1:Xsamp/2))));
X_highMax = find(abs(X(Xsamp/2+1:Xsamp)) == max(abs(X(Xsamp/2+1:Xsamp))))+Xsamp/2;
domsamples = zeros(length(X),1);
domsamples(X_lowMax) = 1;
domsamples(X_highMax) = 1;

domsamples = domsamples.*X;

domsamplestime = ifft(domsamples);
Edom1 = sum(abs(domsamplestime).^2);



% Etot frequency domain
Etot2 = 1/N*sum(abs(X).^2);

Edom2 = 1/N*sum(abs(domsamples).^2);

%%
sound(abs(domsamplestime),8000); % make sure that the quality of the recording is okay
%% Harmonic 

harmonic1 = 1-Edom1/Etot1;
harmonic2 = 1-Edom2/Etot2;

%% AR


mo = ar(x,2);
A = mo.a;
[R,P,K] = residue(1,A);

placement = norm(P(1));
distance = 1-placement;

%% Vowel

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

edata = iddata(aa(1:2*floor(N/3)),[],T);
vdata = iddata(aa(2*floor(N/3)+1:N),[],T);
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

%adf = etfe(aa)

n = length(mo.a);


w=0:1/(1000-1):1;
z = exp(-i.*w);
z_vec = ones(length(z),1);
for i=1:n-1
    for j=1:length(z_vec)
        z_vec(:,1)=z_vec(:,1).*z';
    end
z_vec = [ones(length(z),1) z_vec];
end
 
z_vec_done = z_vec*mo.a';
















