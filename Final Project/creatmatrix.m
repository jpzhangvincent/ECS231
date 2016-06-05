function [A, X, Y] = creatmatrix(n,mu,r,c)
% e.g creatmatrix(8,2,10,8) ouputs 10,8 matrix
X = randn(n,r);
Y = mu+randn(n,c);
A = [];
for i = 1:r
    for j = 1:c
        A(i,j) = log(norm(X(:,i) - Y(:,j)));
    end
end
        

