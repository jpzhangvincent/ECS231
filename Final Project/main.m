% First experiment:
% create a 4000,4000 matrix with mu=1 such that the ratio for the first
% two eigenvalues is large
%[A, X, Y] = creatmatrix(100,1,4000,100);

%k=3

% compute the exact truncated SVD
%[Uk,Sk,Vk] = svds(A,k);
%Ak = Uk*Sk*Vk';
% compute the approximate rank-k SVD
%[C,U,R]= ASCUR(A, k , 1);
%Acur = C*U*R;
% use Error ratio defined as follows to access the accuracy for our
% algorithms
%ErrorRatio = norm(A - Acur,'fro')/norm(A - Ak, 'fro');

%For Randomized SVD algorithms
%with mu = 1 for large eigenvalues ratio; and mu = 2.5 for small eigenvalue
%ratio
%Case 1:
%For fixed k=1, plot the error ratios convergence varying on p, randSVD_SI
%with LU reorthogonalization or QR reorthognalization
% k = 3;
% l = 15; %we tested on l=5,10,15
% [Uk,Sk,Vk] = svds(A,k);
% Ak = Uk*Sk*Vk';
% t1 = [];
% t2 = [];
% eps_ls1 = [];
% eps_ls2 = [];
% for q = 1:20
%      tStart1=tic;
%      [U,S,V]= randSVD_SI(A,k,l,q,false);
%      tElapsed1=toc(tStart1);
%      t1 = [t1 tElapsed1];
%      A_approx1 = U*S*V';
%      tStart2=tic;
%      [U,S,V]= randSVD_SI(A,k,l,q,true);
%      tElapsed2=toc(tStart2);
%      t2 = [t2 tElapsed2];
%      A_approx2 = U*S*V';
%      ErrorRatio1 = norm(A - A_approx1)/norm(A - Ak);
%      ErrorRatio2 = norm(A - A_approx2)/norm(A - Ak);
%      eps_ls1 = [eps_ls1 ErrorRatio1];
%      eps_ls2 = [eps_ls2 ErrorRatio2];
% end
% 
% plot(1:20,eps_ls1,1:20,eps_ls2)
% legend('with QR factorization','with LU factorization')
     
%calculate average running time
% mean(t1) %for using QR factorization  0.0136, 0.0221, 0.0325 s|0.0131 0.0229 0.0355
% mean(t2) %for using LU factorization  0.0127, 0.0212, 0.0255 s|0.0142 0.0204 0.0303

%Case 2: 
% a)
[A, X, Y] = creatmatrix2(4000,4000,100);
k= 10;
tStart=tic;
[Uk,Sk,Vk] = svds(A,k);
tElapsed=toc(tStart1);
Ak = Uk*Sk*Vk';

%b)
[A, X, Y] = creatmatrix(100,1,4000,100);
k= 10;
tStart=tic;
[Uk,Sk,Vk] = svds(A,k);
tElapsed=toc(tStart1);
Ak = Uk*Sk*Vk';

%c)
[A, X, Y] = creatmatrix(100,4,4000,100);
k= 10;
tStart=tic;
[Uk,Sk,Vk] = svds(A,k);
tElapsed=toc(tStart1);
Ak = Uk*Sk*Vk';
%For fixed k=10, plot the error ratios convergence varying on iteration q, comparing 
%randSVD_SI(l=15), randSVD_impSI(l1=20,l2=15,useLU=false), and randomSVD_impOS(p=5,c=5)
eps_ls1 = [];
eps_ls2 = [];
eps_ls3 = [];
t1 = [];
t2 = [];
t3 = [];
l1 = 20;
l2 =15;
p = 5;
c = 5;
tStart=tic;
[Uk,Sk,Vk] = svds(A,k);
tElapsed=toc(tStart1);
Ak = Uk*Sk*Vk';
for q = 0:40
     tStart1=tic;
     [U1,S1,V1]= randSVD_SI(A,k,l,q,true);
     tElapsed1=toc(tStart1);
     t1 = [t1 tElapsed1];
     tStart2=tic;
     [U2,S2,V2]= randSVD_impSI(A,k,l1,l2,q,true);
     tElapsed2=toc(tStart2);
     t2 = [t2 tElapsed2];
     tStart3 = tic;
     [U3,S3,V3]= randSVD_impOS(A,k,q,p,c,true);
     tElapsed3=toc(tStart3);
     t3 = [t3 tElapsed3];
     
     A1_approx = U1*S1*V1';
     A2_approx = U2*S2*V2';
     A3_approx = U3*S3*V3';
     ErrorRatio1 = norm(A - A1_approx)/norm(A - Ak);
     ErrorRatio2 = norm(A - A2_approx)/norm(A - Ak);
     ErrorRatio3 = norm(A - A3_approx)/norm(A - Ak);
     eps_ls1 = [eps_ls1 ErrorRatio1];
     eps_ls2 = [eps_ls2 ErrorRatio2];
     eps_ls3 = [eps_ls3 ErrorRatio3];
end

plot(1:40,eps_ls1(2:41),1:40,eps_ls2(2:41),1:40,eps_ls3(2:41))
legend('Randomized Subspace Iteration','Improved Randomized Subspace Iteration','Improved Randomized with Oversampling')
xlabel('p') % x-axis label
ylabel('Error Ratio for Accuracy') % y-axis label

% 
% %calculate average running time  svds: 187.9236 | 224.1788 | 365.8224
mean(t1) %for using QR factorization  0.0421 | 0.0401 | 0.0442
mean(t2) %for using LU factorization  0.0446 | 0.0435 | 0.0473
mean(t3) %for using LU factorization  0.0485 | 0.0471 | 0.0509


%Case 3:
%For randomSVD_impOS method, plot the error ratios convergence varying on c
%to see how the accuracy changes
%a) Fixed p and q, different c
A = creatmatrix(100,2,4000,4000);
k = 50;
[Uk,Sk,Vk] = svds(A,k);
Ak = Uk*Sk*Vk';
p = 5;
q = 5;
t1 = [];
eps_ls1 = [];

for c = 1:20
     tStart1=tic;
     [U,S,V]= randSVD_impOS(A,k,q,p,c,true);
     tElapsed1=toc(tStart1);
     t1 = [t1 tElapsed1];
     A_approx1 = U*S*V';
     
     ErrorRatio1 = norm(A - A_approx1)/norm(A - Ak);
     
     eps_ls1 = [eps_ls1 ErrorRatio1];
end

plot(1:20,eps_ls1)
%legend('with QR factorization','with LU factorization')
xlabel('c') % x-axis label
ylabel('Error Ratio for Accuracy') % y-axis label

%b) Fixed q and c, different p
c = 5;
q = 5;
t1 = [];
eps_ls1 = [];

for p = 1:20
     tStart1=tic;
     [U,S,V]= randSVD_impOS(A,k,q,p,c,true);
     tElapsed1=toc(tStart1);
     t1 = [t1 tElapsed1];
     A_approx1 = U*S*V';
     
     ErrorRatio1 = norm(A - A_approx1)/norm(A - Ak);
     
     eps_ls1 = [eps_ls1 ErrorRatio1];
end

plot(1:20,eps_ls1)
%legend('with QR factorization','with LU factorization')
xlabel('p') % x-axis label
ylabel('Error Ratio for Accuracy') % y-axis label

%Case 4 real application   work
A = Problem.A;
A = full(A);
k =10;
tStart=tic;
[Uk,Sk,Vk] = svds(A,k);
tElapsed=toc(tStart);
Ak = Uk*Sk*Vk';
eps_ls1 = [];
eps_ls2 = [];
t1 = [];
t2 = [];
l1 = k+20;
l2 =k+10;
p = 10;
c = 10;
for q = 0:20
     tStart1=tic;
     [U1,S1,V1]= randSVD_impSI(A,k,l1,l2,q,true);
     tElapsed1=toc(tStart1);
     t1 = [t1 tElapsed1];
     tStart2 = tic;
     [U2,S2,V2]= randSVD_impOS(A,k,q,p,c,true);
     tElapsed2=toc(tStart2);
     t2 = [t2 tElapsed2];
     
     A1_approx = U1*S1*V1';
     A2_approx = U2*S2*V2';
     
     ErrorRatio1 = norm(A - A1_approx)/norm(A - Ak);
     ErrorRatio2 = norm(A - A2_approx)/norm(A - Ak);
     
     eps_ls1 = [eps_ls1 ErrorRatio1];
     eps_ls2 = [eps_ls2 ErrorRatio2];
    
end

plot(1:20,eps_ls1(2:21),1:20,eps_ls2(2:21))
legend('Improved Randomized Subspace Iteration','Improved Randomized with Oversampling')
xlabel('p') % x-axis label
ylabel('Error Ratio for Accuracy') % y-axis label
