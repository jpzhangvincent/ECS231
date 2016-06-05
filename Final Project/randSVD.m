function [U, S, V]= randSVD(A,k,l)
[m, n ] = size(A);
if m < n,
    disp('The input matrix should have number of rows m larger than or equal to number of column n');
    return
end
if l < k,
    disp('l should be larger than k');
    return
end 
O = randn(n, l);
Y = A*O;
Q = orth(Y);
B = Q'*A;
[W,S,V] = svd(B,'econ');
U = Q*W;

U = U(:,1:k);
S = S(1:k,1:k);
V=V(:,1:k);

%Ak = U*S*V;