function [J] = jacobiana( fname,x )
%Matriz Jacobiana de fname: R^n --> R^m
%   matriz J de orden mxn
%--------------------------------------
n=length(x);
r=feval(fname,x);
m= length(r);
J=zeros(m,n);
h=1.e-05;

for k=1:n
    xp=x; xp(k)=xp(k) +h;
    J(:,k)=(feval(fname,xp)-r)/h;
    
end


end

