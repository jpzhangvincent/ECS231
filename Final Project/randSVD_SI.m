function [U, S, V]= randSVD_SI(A,k,l,q,useLU)
%inputs: A is a m,n matrix
%        k the target rank
%        l the interger to specify the number of column for a random
%         start matrix , l>=k
%        q number of iteration for power methods
[m, n] = size(A);
if m < n,
    disp('The input matrix should have number of rows m larger than or equal to number of column n');
    return
end
if l < k,
    disp('l should be larger than k');
    return
end 
X = randn(n, l);
Y = A*X;
if(~useLU)
    [Q,R] = qr(Y, 0 ) ;
    for j = 1 : q
        Y = A'*Q;
        [Q,R] = qr(Y, 0 ) ;
        Y = A*Q;
        [Q,R] = qr(Y, 0 ) ;
    end
else 
    [Q,R] = lu(Y);
    for it = 1:q
      Q = (Q'*A)';
      [Q,R] = lu(Q);
      Q = A*Q;
      if(it < q)
        [Q,R] = lu(Q);
      end
      if(it == q)
        [Q,R] = qr(Q,0);
      end
    end
end
clear R;

B = Q'*A;
[W, S ,V] = svd(B,'econ');
U = Q*W;

U = U(:,1:k);
S = S(1:k,1:k);
V=V(:,1:k);
%Ak = U*S*V';