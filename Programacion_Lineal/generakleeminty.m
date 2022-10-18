function [A,b,c] = generakleeminty(n)
c = ones(1,n);
c = -1*c;
b= zeros(n,1);
b(1) = 1;
for i=2:n
    b(i) = 2.^i-1;
end
    A = tril(2*ones(n,n),0) - eye(n);
    A(1,1) = 1;
end