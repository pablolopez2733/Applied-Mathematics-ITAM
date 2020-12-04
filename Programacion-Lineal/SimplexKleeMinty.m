% Prueba Kleeminty
for i=2:10
    tic
    [A,b,c] = generakleeminty(i);
    b=b';
    [x0,z0,ban,numiter] = mSimplexFaseII(A,b,c);
    numiter;
    
    char = 'Para m= %1.0f se tarda %4.4f cpu time y termina en %4.0f interaciones';
    t= toc;
    fprintf(char,i,t,numiter)
    fprintf('\n')
end