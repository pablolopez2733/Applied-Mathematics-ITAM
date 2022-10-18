% TAREA 2
%   Verifiquen su codigo cin un problema que hayan resuelto
%   Este problema lo resolvimos en clase
%   Solucion de clase x=(1,2,0)    valor optimo = 11

A  = [  1    2    3;
        2    2    1 ];
    
b  = [ 5; 6];
c  = [ 3; 4; 5]; % costos
[xo, zo, ban, iter]=mSimplexDualtabl(A,b, c)

% Problema no acotado

A  = [  -1    1/2    -1;
         1     -1    -1 ];
    
b  = [ 3; 2];
c  = [ 2; 1; 2]; % costos
[xo, zo, ban, iter]=mSimplexDualtabl(A,b, c)
