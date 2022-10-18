function [gx, outgx] = restricciones(x)
% 
%
n = length(x)/2;

X = [x(1:n);pi;x(n+1:2*n);0];

n = length(X)/2;
%disp(n);

gnon = zeros(0.5*(n*(n-1)),1);
glin = zeros(n-1,1);
gbounds = zeros(4*n,1);

th = X(1:n); r = X(n+1:2*n);

for i=1:n-1
    for j=i+1:n
        gnon(i*(n-1)+j) = 1 - (r(i)^2 + r(j)^2 - 2*r(i)*r(j)*cos(th(i) - th(j)))-.2;
    end
end

for i=1:n-1
    glin(i) = th(i+1) - th(i);
end

gbounds(1:2*n) = X;
gbounds(2*n+1:3*n) = pi - th;
gbounds(3*n+1:4*n) = 1 - r;
gbounds(4*n+1:5*n) = r;


gx = [gnon
    glin
    gbounds];

if nargout > 1
        outgx = struct(...
            'gnon', gnon,...
            'glin', glin,...
            'gbounds', gbounds);
end
        