function [x, fx, output] = puntosInteriores(...
    fobj,fineq, x0, SIGM)
%              METODO DE PUNTOS INTERIORES
%
% Se programa el metodo de puntos interiores para el problema de
% minimizacion:
%
%               min:    fobj(x)
%               s.t  ineqc(x) >= 0
%
%
% utilizacion:
%          [x, fx, output] = puntosInteriores(fobj, ineqc, x0)
%
% variables de entrada:
%              fobj := funcion objetivo (string)
%             ineqc := funcion de restricciones (string)
%                x0 := punto inicial
%              gam0 := (opcional) valor inicial de gamma.
%
% variables de salida:
%                x := punto donde se encontro el minimo
%               fx := (opcional) valor de la funcion objetivo en x.
%           output := (opcional) estructura que muestra varios
%                     parametros:
%                       - iteraciones
%                       - stepsize
%                       - CNPO
%                       - constviolation
%                       - mLagrange
%                       - trayectoria
%
if nargin < 4
    SIGMA = 0.3;
else
    SIGMA = SIGM;
end

gamk = 1;
tol = 1e-6;
MaxIter = 75;
p = length(feval(fineq,x0));
n = length(x0);

xk = x0;
e = ones(p,1);
idx = 1:p;

zk = e;
muk = e;
Zk = eye(p);
Uk = Zk;
invZk = diag(1./zk);
invUk = diag(1./muk);

gk = feval(fineq,xk);
Ck = jacobiana(fineq,xk);
gradfk = gradiente(fobj,xk);

Bk = eye(n);

gradLk = -gradfk + Ck'*muk;
epsik = -zk + gk;
sigmak = -muk.*zk + gamk;

Trayectoria = zeros(n,MaxIter);

termine = false;
k=0; kfail =0; ii = 0;

while ~termine 
   
    deltaXk = linsolve(Bk + Ck'*invZk*Uk*Ck , ...
       gradLk - Ck'*invZk*Uk*epsik + Ck'*invZk*sigmak);
    
    deltaMuk = - invZk*Uk*(epsik - invUk*sigmak + Ck*deltaXk);
    deltaZk = invUk*(sigmak - Zk*deltaMuk);
    
    idMu = idx(deltaMuk<0);
    idZ = idx(deltaZk<0);
    
    if ~isempty(idMu)
        alphaMuk = min(-muk(idMu)./deltaMuk(idMu));
    else
        alphaMuk = 1;
    end
    
    if ~isempty(idZ)
        alphaZk = min(-zk(idZ)./deltaZk(idZ));
    else
        alphaZk = 1;
    end
    
    alphak = 0.995*min([1 alphaMuk alphaZk]);
    
    xk = xk + alphak*deltaXk;
    muk = muk + alphak*deltaMuk;
    zk = zk + alphak*deltaZk;
    
    yk = -gradfk; % preliminar
    DiffHk = - Ck; % preliminar
    
    if ~isempty(zk(zk==0))
        zk(zk==0) = 1e-6;
    end
    if ~isempty(muk(muk==0))
        muk(muk==0) = 1e-6;
    end
    
    Zk = diag(zk);
    Uk = diag(muk);
    
    invZk = diag(1./zk);
    invUk = diag(1./muk);
    
    gk = feval(fineq,xk);
    Ck = jacobiana(fineq,xk);
    gradfk = gradiente(fobj,xk);
    
    DiffHk = Ck + DiffHk;
    
    sk = alphak*deltaXk;
    yk = yk + gradfk + DiffHk'*muk;
    
    if cond(Bk)>10000
        Bk = eye(n);
        kfail = kfail + 1;
    end
    
    skyk = sk'*yk; 
    
    if skyk == 0
        if norm(sk) < 1e-16
            ii = ii + 1;
        end
        %disp(norm(deltaXk));
        %disp(k);
        skyk = 1e-4;
    end
    
    Bs = Bk*sk;
    sBs = sk'*Bs;
    %SIGMA = 0.26;
    
    if skyk <= SIGMA*sBs ;
        THk = (1-SIGMA)*(sBs/(sBs - skyk));
    else
        THk = 1;
    end
    yk = THk*yk + (1-THk)*Bs;
    

    Bk = Bk - (Bs*Bs')/sBs + (yk*yk')/skyk;
    
    gamk = gamk/10;
%    tauk = tauk*0.5;
    Trayectoria(:,k+1) = xk;

    gradLk = -gradfk + Ck'*muk;
    epsik = -zk + gk;
    sigmak = -muk.*zk + gamk;
    
    k = k+1;
    Fk = [gradLk;epsik;muk.*zk];
    termine = norm(Fk)<=tol || ...
        k>=MaxIter;
    %|| ii >= 30  
    
end

Trayectoria = [x0 Trayectoria];
Trayectoria(:,k:MaxIter) = [];

x = xk;
if nargout > 1
    fx = feval(fobj,x);
    if nargout > 2
        output = struct(...
            'iteracionesTotales', k,...
            'kfail',kfail,...
            'stepsize', norm(alphak*[xk;muk;zk]),...
            'CNPO', norm(Fk), ...
            'F', Fk, ...
            'z', zk,...
            'mLagrange',muk,...
            'trayectoria', Trayectoria);
    end
end
