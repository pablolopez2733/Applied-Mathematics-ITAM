function [d] = Perimetro(x)
    m = length(x);
    n = m/2;
    theta = x(1:n);
    r = x(n+1:m);
    d = r(1) + r(n);
    for k=1:n-1
        h=sqrt(r(k)^2+r(k+1)^2-2*r(k)*r(k+1)*sin(theta(k)-theta(k+1)));
        d=d+h;
    end
    d=-d;
end