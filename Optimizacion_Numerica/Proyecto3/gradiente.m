function [g]= gradiente(fname,a)
%
h = 1.e-05;
fa = feval(fname,a);
n = length(a);
g = zeros(n,1);
at=a;
for i=1:n
    at(i) = at(i) + h;
    fat = feval(fname,at);
    g(i)= (fat-fa)/h;
    at(i)=a(i);
end