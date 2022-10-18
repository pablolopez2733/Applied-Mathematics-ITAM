function x0 = calculatex0(original,A,B)

[m,n] = size(original);

check = (1:n);

x0 = (1:n);

for i = 1:n
	if (ismember(i,B) && ismember(i,check))
		for j = 1:n
			if(B(j) == i)
				x0(i) = A(j,end);
			end
		end
	else
		x0(i) = 0;
	end
end

end
