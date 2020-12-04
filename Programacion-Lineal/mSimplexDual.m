function [xo, zo, ban, iter] = mSimplexDual(A, b, c)
% Marcos Juan de Dios
% In : A ... mxn matrix
%purpose: Versión del Simplex en la Fase II minimizar c^T x
%sujeto a Ax >= b , x >= 0 , c >= 0
%       b ... column vector with as many rows as A
%       c ... row vector with as many columns as A
% xo ... SFB  optima del problema zo ... valor  ?optimo del problema ban ... indica casos:
%      -1 ... si el conjunto factible es vacio
%       0 ... si se encontro una solucion  optima
%       1 ... si la funcion objectivo no es acotada.
% iter ... es el numer ?o de iteraciones (cambios de variables basicas) que
% hizo el método

[m,n]=size(A);
xn=zeros(n,1);
x=[xn;-b];
A=[-A, eye(m)];
c=[c;zeros(m,1)];
% cn=c;
% cb=zeros(m,1);

b=-b;
h=b;
N=1:n;
B=n+1:n+m;
iter=0;
ban=0;                          % Suponemos que el problema tiene sol

while min(h)<0
    H=A(:, B)\A(:, N);          % Equivale a usar inv(A(:, B))*A(:, N)
    lambda= (A(:, B)')\(c(B));  % Equivale a usar ((c(B)')*inv(A(:, B)))'
    z=lambda'*A;                
    
    ind=(find(h'<0));
    sale=min(B(ind));           % Encontramos la variable de Entrada
    s=find(B==sale);
    
    if 0 <= min(H(s,:))         % Condicion necesaria para que el Dual no 
        ban = -1;               % sea acotado, implica que Cf(P) es vacio
        break;
    end
    vecPiv=(z(N)-c(N')')./H(s,:);   
    ind=(find(H(s,:)<0));
    piv=min(vecPiv(ind));       % Criterio para escojer var de entrada
    e=find(vecPiv==piv);        % Lugar donde esta la var de entrada
    entra=N(e);                 % Encontramos la variable de entrada
    
    B(B==sale) = entra;         % Actualizamos variable Basicas
    N(N==entra) = sale;         % Actualizamos variable No-Basica
    
    h=A(:, B)\b;                % Equivale a usar inv(A(:, B))*b
    iter = iter + 1;
end
x=zeros(n+m,1);
x(B)=h;
xo=x(1:n);
zo=c'*x;
    
% PROBLMEA 1
%   Ban = 1 nunca puede ocurrir, ie el problema no va a estar no acotado.
%   Pues por la construcion del problema primal Ax >= b , x >= 0 , c >= 0
%   el problema dual es A'y <= c , y >=0 y c >=0, y como c >= 0 entonces 
%   el vector de ceros 0 va a pertenecer a Cf(D). El teorema dice que si el
%   problmea P no esta acotado, entonces el  Cf(D) es vacio, por lo tanto
%   podemos afirmar que el problema P no esta no acotado (ie
%   ban~=1 para los problemas que estamos tratando)
        


