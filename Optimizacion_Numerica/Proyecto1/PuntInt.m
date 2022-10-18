%Proyecto Optimización Numérica  
%Eric Bazaldua Miñana
%Francisco Gabriel Huerta Fernández 166040
%Pablo López Landeros 
% Codigo 1: Puntos interiores para resolver problema de programación
% cuadrática 
% Resolver min .5x^T*Q*x + c^T*x
%           s.a Ax>=b 
% Con Q spd y de (n*n)
%     c vector de (n*1)
%     A matriz de (m*n) (matriz de restricciones)
%     b vector de (m*1) (vector de restricciones)
% Función PuntInt; 
%entrada: Q,A,c,b x0,y0,m0,sigma
% x0 es el valor inicial tal que Ax0>b
% y0 y m0 son vectores iniciales tal que y0,m0>0
% sigma en [0,1] 
%Sálida: 
%vector x: aproximación del mínimo local
%vector y: variable de holgura en la restricción 
%vector m: aproximación del multiplicador de lagrange asociado a la
%restricción. 
function[x,y,m]= PuntInt(Q,A,c,b)
%Primero hicimos un paso 
sigma=.5;
m=length(b); 
n=length(c); 
x=ones(n,1); %fijamos x0 empezamos en (1,....,1)^T
y=ones(m,1); %fijamos y0 
mu=ones(m,1); %fijamos m0, ambos satisfacen que y0,m0>0
U=diag(mu);
Y=diag(y);
e=ones(m,1);
na1=y'*mu/m; %el valor nabla'
na=sigma*na1; %valor nabla
F1=Q*x-A'*mu+c; 
F2=A*x-y-b;
F3=Y*U*e-na*e;
F=[F1;F2;F3]; %función F(x,y,mu)
tol=1e-8;
cont=0; %contador 
max=100;
while (norm(F)>0.00000001 && cont<max) 
    %realizamos la iteraciones necesarias mientras se satisfaga que
    %||F(x,y,mu)||_2< tol para cierta tolerancia tol
    U=diag(mu);
    Y=diag(y);
    e=ones(m,1);
    na1=y'*mu/m;
    na=sigma*na1;
    rx=c-A'*mu+Q*x;
    ry=A*x-y-b;
    rm=Y*U*e-na*e;
    P=(Q+A'*inv(Y)*U*A);
    %P1=rx+A'*inv(Y)*U*(ry+y-na*inv(U)*e);
    %Pchol=chol(P,'lower') %usamos cholesky para obtener el valor de Dx porque la expresión P es complicada 
    %Dx=Pchol\(Pchol'\P1);
    %[Dx,~,~,~,~] = pcg(P,P1,tol);
    Dx=-inv(P)*(rx+A'*inv(Y)*U*(ry+y-na*inv(U)*e));
    Dy=A*Dx+ry;
    Dm=-inv(Y)*(U*Dy+rm);
    k1=length(Dy);
    k2=length(Dm);
    
    for i= 1:k1
        if Dy(i)>=0
            g1(i)=1;
        else
            g1(i)=-y(i)*(1/Dy(i));
        g1;    
        end
    end
     for i= 1:k2
        if Dm(i)>=0
            g2(i)=1;
        else
            g2(i)=-mu(i)*(1/Dm(i));
        g2;    
        end    
     end
    a1=min(g1);
    a2=min(g2);
    a3=min(a1,a2);
    a=.995*min(1,a3);
    %loop 
    x=x+a*Dx; %valor de regreso para aproximación de x
    y=y+a*Dy; %valor de regreso para aproximación de y
    mu=mu+a*Dm; %valor de regreso para aproximación de mu
    F1=Q*x-A'*mu+c;
    F2=A*x-y-b;
    F3=Y*U*e-na*e;
    F=[F1;F2;F3];
    cont=cont+1; %conteo final de pasos
    if mod(cont, 2) == 0
        sigma = 0.0;            % Paso de predicción
    else
        sigma = 0.5;            % Paso de correción
    end
end
end 