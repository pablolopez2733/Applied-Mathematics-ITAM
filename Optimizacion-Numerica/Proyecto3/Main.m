% Eric Bazaldua Miñana 155279
% Pablo Lopez Landeros 178863
% Francisco Gabriel Huerta Fernandez 166040
% Script para calcular perimetro máximo de poligonos. Se necesita definir
% la cantidad de vertices deseadas
clear all
close all
warning('off','all');
% ----- ELIGE EL NUMERO DE VERTICES --------- % 
n = 10;
x = [];
y = [];
x_matlab = [];
y_matlab = [];
%Creamos el vector inicial
x0 = puntoInicial(n);

%Solución por el método puntos interiores
[xStar,muStar,iter] = puntosInteriores('Perimetro','restricciones',x0);
% Solucion matlab
restPoligono = @(x) deal(-restricciones(x),[]);
options = optimset('Algorithm','interior-point','Display','off');
[xStarMatlab,fx,exitflag,output] = fmincon('Perimetro',x0,[],[],[],[],[],[], restPoligono,options);

%Puntos necesarios para la gráfica
[x,y]=pol2cart(xStar(1:length(xStar)/2),xStar(length(xStar)/2+1:length(xStar)));
[x_matlab,y_matlab]=pol2cart(xStarMatlab(1:length(xStarMatlab)/2),xStarMatlab(length(xStarMatlab)/2+1:length(xStarMatlab)));
x=[0;x];
y=[0;y];
x_matlab=[0;x_matlab];
y_matlab=[0;y_matlab];

%Grafica de nuestra Solución vs Solución matlab
plot([x; x(1)],[y;y(1)],'b');
hold on;
plot([x_matlab; x_matlab(1)],[y_matlab;y_matlab(1)],'r');
legend('Solución nuestra','Solución Matlab')
lados = round(n,0);
title(sprintf('Polígono con %g lados', lados))