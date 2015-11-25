function [xhat, w2] = decdown(xn, p, N)

xhat = xn(1:p:length(xn));

w2 = (0:N/p-1)*2*pi/(N/p);

end

