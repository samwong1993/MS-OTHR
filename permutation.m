function PMat = permutation(M)
A = rand(M,1);
B = A(randperm(M));
[~,IA] = sort(A);
[~,IB] = sort(B);
I(IB) = IA;
PMat(:,I) = eye(length(A));
end




