function [U, S, V]= randSVD_impSI(A,k,l1,l2,q,useLU)
%inputs: A is a m,n matrix
%        k the target rank
%        l1 the interger to specify the number of column for the first random
%          start matrix , l1>l2>=k
%        l2 the interger to specify the number of column for the second random
%          start matrix , l1>l2>=k
%        q  number of iteration for power methods
%        useLU boolean value to specify to use LU factorization to do
%        reorthogonalization

[m, n] = size(A);
if m < n,
    disp('The input matrix should have number of rows m larger than or equal to number of column n');
    return
end
if l2 < k||l1<l2
    disp('l2 should be larger than k and l1 should be larger than l2');
    return
end
[U1,S1,V1] = randSVD(A,l2,l1);
[U, S, V] = truncatedSVD_SI(A,k,l2,q,V1,useLU);

%Ak = U*S*V;