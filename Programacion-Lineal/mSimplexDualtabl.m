function [x0,z0,ban,iter] = mSimplexDualtabl(A, b, c)

[m1 , n1] = size(A);
original = A; % Gardamos original para calculos de x0, funcion externa
A = [(-1)*A,(1)*eye(m1)]; % Multiplicar A por -
aux = zeros(m1);
tableau = [A;c aux(1,:)]; % Queda el vector c encima de la matriz A
b = [b,0]; % Agregando un 0 al primer valor de c para el tableau
b = (-b)'; % Multiplicar b por -1
tableau = [tableau b]; % Finalmente, agregamos b a tableau
[m2 , n2] = size(tableau); % Calculamos tamano del tableau, para uso despues
%[A,p,q,e,B] = ChoosePivotDual(tableau,tableau(3:4)); % Dual en general.
% Importantemente, regresa e, el valor del pivote, y (p,q), su posicion.
A = tableau;
B = [(m1+1):(m1+n1)];
auxban = 0;
iter = -1;

while auxban == 0
	[A,p,q,e,B,auxban,ban] = ChoosePivotDual(A,B);
	A = simplex(A,p,q);
	iter = iter + 1;
end

z0 = (-1)*A(end,end);

x0 = calculatex0(original,A,B);
