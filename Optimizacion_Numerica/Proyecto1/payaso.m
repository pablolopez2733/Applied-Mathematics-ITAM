ks = [5, 20, 30, 60, 80];
%ks = [1,2,3,4,5];
t_qp = zeros(5,1);
t_pint = zeros(5,1);
for i=1:length(ks)
    tic
    [~,~]=descenso2pasos_qp(X,ks(i));
    t_qp(i)= toc;
    
    tic
    [~,~]=descenso2pasos(X,ks(i));
    t_pint(i)=toc;
end
