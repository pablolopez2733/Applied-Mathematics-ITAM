function [xo, zo, ban, iter] = mSimplexFaseII(A, b, c)
% Comienza agregando variables de holgura:
A=[A,eye(length(b))];

% A continuacion se definen los elementos y las matrices que se utilizaran:
[m,n] = size(A);
B = eye(m); 
cN = [c,zeros(1,length(b));]; % Se agregan ceros en las entradas basicas del vector de costos
cB(1:m) = 0;
l=0;
iter = 0;

% Checa si conjunto factible es vacio
if (any(b<0))
    zo = nan;
    xo = nan;
    ban = -1;
    iter = 0;
    return
end

% Determina las variables básicas y no básicas
for i=1:n
    XN(i) = i;
end
for i = n+1:m+n
    XB(i-n) = i;
end

% Comienza el loop

while(l~=1)
    
% Determina la variable de entrada
lambda = transpose(B) \ transpose(cB);
p=transpose(lambda)*A;
rN=transpose(p) - transpose(cN);

i=1;
while(rN(i)<=0 && i<n)
      i = i+1;
end

if(rN(i)<=0) % Ya se llego al optimo
      zo = transpose(lambda)*transpose(b);
      x = B\transpose(b);
      for i=1:n
            for j=1:length(XB)
                if (XB(j) == i)
                    xo(i) = x(j);
                end
            end
        end
      ban=0;
      return
else
      iter=iter+1;
      xe=XN(i); 
      e=i;
end

% Determina la variable de salida
h=B\transpose(b);
H = B\A(1:m,e);
for i=1:m
      if (H(i) <=0 )
            aux(i) = 1000000; % No participa
      else
            aux(i)=h(i)/H(i);
      end
end
if (min(aux) == 1000000) % Problema no acotado
      xo=nan;
      zo = nan;
      ban = 1;
      return
else
      i =1;
      while(min(aux) ~= aux(i) &&i<=length(aux))
            i = i+1;
      end
      s=i;
      xs=XB(i);
end

% Calcula B+ y N+
A_1=A(1:m,e);
B_1=B(1:m,s);
A(1:m,e)=B_1;
B(1:m,s)=A_1;
hip1=cB(s);
hip2=cN(e);
cB(s)=hip2;
cN(e)=hip1;

% Actualiza las variables basicas y no basicas
i = 1;
    while (XB(i) ~=xs && i<length(XB))
        i = i + 1;
    end
    XB(i) = xe;
    e = i;
    i = 1;
    while (XN(i) ~= xe && i<length(XN))
        i = i + 1;
    end
    XN(i) = xs;
    s=i;
end
end


