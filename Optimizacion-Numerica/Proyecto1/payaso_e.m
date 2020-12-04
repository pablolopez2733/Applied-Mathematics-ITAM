%% Eric Bazaldua Miñana 155279
%% Mauricio Huerta
%% Pablo Torre
% Script que corre utiliza los codigos descenso2pasos.m, descenso2pasos_qp
% y PuntInt para genera la imagen clown
% PuntInt.m: Código para aproximar por puntos interiores 
% descenso2pasos.m: Código de descenso en 2 pasos utilizando puntos
% interiores
% descenso2pasos_qp.m: Código de descenso en 2 pasos utilizando quadprog

ks = [5, 20, 30, 60, 80]; %Vector de las iteraciones 
%ks=[1,2,3,4,5];
t_qp = zeros(5,1); %Vector de los tiempos del quadprog
t_pint = zeros(5,1); %Vector de los tiempos de puntos interiores
load clown
for i=1:length(ks)
    tic
    [W,H]=descenso2pasos_qp(X,ks(i)); %Corremos descenso 2 pasos con quadprog
    t_qp(i)= toc; %Guardamos los tiempos
    fprintf('\n\t%s \t%d','Imagen del payaso utilizando quadprog para k=',ks(i))
    colormap(gray)
    image(W*H)
    pause
    
    tic
    [W,H]=descenso2pasos(X,ks(i)); %Corremos descenso 2 pasos con puntos interiores
    t_pint(i)=toc; %Guardamos los tiempos 
    fprintf('\n\t%s \t%d','Imagen del payaso utilizando puntos interiores para k=',ks(i))
    colormap(gray)
    image(W*H)
    pause
end


%Genera la Tabla de los tiempos de punto interior y quadprog
fprintf('Genera la Tabla de los tiempos de punto interior y quadprog')
fprintf('\n\t%s \t\t%s \t\t\t%s\n' ,'k','quadprog','Puntos Interiores');
fprintf('\t-------------------------------------------------------------------------------------------------');

fprintf('\n\t%s \t\t%d \t\t%.4d', '5',t_qp(1),t_pint(1));
fprintf('\n\t%s \t\t%d \t\t%.4d', '20',t_qp(2),t_pint(2));
fprintf('\n\t%s \t\t%d \t\t%.4d', '30',t_qp(3),t_pint(3));
fprintf('\n\t%s \t\t%d \t\t%.4d', '60',t_qp(4),t_pint(4));
fprintf('\n\t%s \t\t%d \t\t%.4d', '80',t_qp(5),t_pint(5));
fprintf('\n');
fprintf('\t-------------------------------------------------------------------------------------------------');
fprintf('\n');

