function [W, H] = descenso2pasos(X,k)
%% Método de descenso en 2 pasos (actualizando H y W con puntos interiores):
%ENTRADAS:
%X -> matriz de mxn
%k -> constante de factorización

%SALIDAS:
%W -> matriz de mxk
%H -> matriz de kxn

%El resultado se grafica con: image(W*H)

%------------------------------------------------------------------------
%% Codigo

%Definimos dimensiones
[m,n] = size(X);
%Las primeras aproximaciones de W y H tienen que ser positivas, rand genera
%números en (0,1), por lo tanto, lo podemos utilizar.
W = rand(m,k);
H = rand(k,n);

%PARÁMETROS:
maxiter = k;
tol = 1e-06;
i=0; %contador

%Podríamos usar una H y W auxiliar para ir almacenando las iteraciones
%anteriores y poder ver cuanto "avanzamos" en cada iteración:
H_aux=zeros(k,n);
W_aux=zeros(m,k);

%inicializamos gen_step que será lo que se comparará con la tolerancia
%gen_step = ||Hk+1 - Hk||_F + ||Wk+1 - Wk||_F
gen_step=10; 


while (i < maxiter && gen_step > tol) 

    %Resolvemos para H:
    Q=W'*W;
    c=W'*X;
    A = eye(k);
    b= zeros(k,1);
    
    for j=1:n
        H_aux(:,j)=H(:,j);% en la aux, guardo la anterior
        [rH,y,na]= PuntInt(Q,A,-c(:,j),b); %resolvemos el problema
        step_H = H_aux(:,j)-rH';%esta es la diferencia en cada paso
        H(:,j)=rH';
    end
    
    %Resolvemos para W:
    Q=H*H';
    c=H*X';
    A = eye(k);
    b= zeros(k,1);
  
    for u=1:m
        W_aux(u,:)=W(u,:);
        [rW,y,ma]= PuntInt(Q,A,-c(:,u),b);
        W(u,:)=rW;
        step_W = W_aux(u,:)-rW;
    end
    
i=i+1;

%cada iteracion calculamos gen_step para comparar con la tol
gen_step = norm(step_H,'fro')^2 + norm(step_W,'fro')^2;
end

%% Para probar:
% load clown
% colormap(gray)
% image(X)

end

