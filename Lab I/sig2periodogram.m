function [ Phi,f ] = sig2periodogram( x,T )
%SIG2PERIODOGRAM s.101 in the signal processing book

if nargin <2 , T=1; end
N=length(x);
X=fft(x);
Phi=X.*conj(X);
Phi=Phi(1:round(N/2)0+1)/N*T;
f=(0:round(N/2))/N/T;
end

