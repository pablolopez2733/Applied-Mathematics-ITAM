function [A,p,q,e,B,auxban,ban] = ChoosePivotDual(A,B)
% input: augmented matrix A, basis matrix B
%
% output : augmented matrix A,
%         pivot column q , pivot row p , pivot element e ,
%         basis matrix B
ban = 0;

[m,n] = size(A);
min=0;
p=0;
auxban = 0;

for I=1:(m-1)
  if A(I,end)< min
    min = A(I,n);
    p=I;
  end
end

if p == 0
  auxban = 1;
  disp('solucion optima encontrada')
  ban = 0;
  e=0;
  q=0;
  return;
end

min = Inf;
q = 0;
count = 0;
for k = 1:n-1
   if A(m,k) ~= 0
     if A(p,k) < 0
       col = abs(A(m,k)/A(p,k));
       if col < min
        min = col;
        q = k;
       end
    end
   end
end

for k = 1:n-1
  if A(p,k) >= 0
    count = count +1;
  end
end

if count == n-1
  auxban = 1;
  disp('problema no acotado')
  ban = 1;
  e=0;
  q=0;
  return
end

e = A(p,q);
B(p) = q;
return
