function C= NOCS(A, k, c)
n = size(A,2);
%c = 2*k/epsi*(1+0.01);
epsi = 2*k/c*(1+0.01);
[U, S, V]= randProjSVD(A, k , 3 );
U = A - U*S*V;
V = V';
s = DSS(V,U,round(c - 2*k/epsi));
C1 = A*diag(s);
C1( :, all( ~any( C1 ), 1 )) = [];
D = A - C1*pinv(C1)*A;
p = [];
normd2 = norm(D)^2;
for i = 1:n
    p(i) = D(:,i)'*D(:,i)/normd2;
end
c2 = 2*k/epsi;
c2 = round(c2);
index = AdaptiveSampling(p, c2);
C2 = A(:,index);
C = [C1,C2];

end
    