% Prueba Kleeminty
for i=2:10
    tic
    [A,b,c] = generakleeminty(i);
    b=b';
    [x0,z0,ban,iteraciones] = mSimplexFaseII(A,b,c);
    iteraciones;
    
    formatSpec = 'Para m= %1.0f se tarda %4.4f cpu time y termina en %4.0f interaciones';
    t= toc;
    fprintf(formatSpec,i,t,iteraciones)
    fprintf('\n')
end