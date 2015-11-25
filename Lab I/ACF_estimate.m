function ACF = ACF_estimate( y )
%Bartleys ACF estimet

N = length(y);
ACF = zeros(1,N);
for k = -floor(N/2):floor(N/2)
    for n = 1:(N-abs(k))
        ACF(k+floor(N/2)+1) = ACF(k+floor(N/2)+1)+y(n)*y(n+abs(k))./N;
    end
end
end