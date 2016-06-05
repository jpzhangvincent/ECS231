function [U, S, V]= randSVD_impSI(A,k,l1,l2,q,useLU)
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