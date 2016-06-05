function [C,U,R]= ASCUR(A, k , epsi )
n = size(A,2);
c = 2*k/epsi;
r = c/epsi*(1+epsi);
C= NOCS(A, k, c);
r1 = c;
R1t = NOCS(A',k, r1);
R1 = R1t';
r2 = c/epsi;
D = A - A*pinv(R1)*R1;
p = [];
normd2 = norm(D)^2;
for i = 1:n
    p(i) = D(:,i)'*D(:,i)/normd2;
end
r2 = round(r2);
index = AdaptiveSampling(p, r2);
R2 = A(index,:);
R = [R1',R2']';
U = pinv(C)*A * pinv(R);


