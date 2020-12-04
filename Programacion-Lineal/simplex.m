function A = simplex(A,p,q)
% Input:	  Augmented matrix A, pivot row p, pivot column q
% Output:   Augmented matrix A

[m,n] = size(A);
format rat

if (p ~= 0 || q ~= 0)
	C = zeros(m,n);
	J = (1:n);

	C(p,J)=(A(p,J)/A(p,q));

	for I = 1:m
		if I ~= p
			C(I,J) = (A(I,J)) - (A(I,q))*C(p,J);
		end
	end

	A = C;

	return
else
	A = A;
end
