function [X0] = puntoInicial(n)
%              GENERA PUNTO INICIAL
%             poligono de area maxima
%
%  
%
t = (1/(2*n):1/n:1)*2*pi - 1/(2*n);

x = 1.*cos(t);
y = 1.*sin(t)+1;

[th,r] = cart2pol(x,y);

dif = abs(th - pi);

[~,k] = min(dif);

if k < n
    th = [th(k+1:n) th(1:k-1) pi];
    r = [r(k+1:n) r(1:k-1) 0];
else
    th = [th(1:n-1) pi];
    r = [r(1:n-1) 0]; 
end

X0 = [th(1:n-1) 0.5.*r(1:n-1)]';