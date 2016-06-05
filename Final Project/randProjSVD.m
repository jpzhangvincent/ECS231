function [U, S, V]= randProjSVD(A, k , q )
[m n ] = size(A);
O = randn( n , k );
Y = A*O;
[Q R] = qr(Y, 0 ) ;
for j = 1 : q
Y = A'*Q;
[Q R] = qr(Y, 0 ) ;
Y = A*Q;
[Q R] = qr(Y, 0 ) ;
end
B = Q'*A;
[U, S ,V] = svd(B);
U = Q*U;
%Ak = U*S*V;