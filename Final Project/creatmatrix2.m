function [A, X, Y] = creatmatrix2(n,r,c)
% e.g creatmatrix2(200,10,10) ouputs 10,10 matrix
X = [];
Y = [];
A = [];

th = 0:pi/n:2*pi;
x1 = sqrt(2) * cos(th) - 1;
y1 = sqrt(2) * sin(th) - 1;
X(1,:) = x1;
X(2,:) = y1;

x2 = 2*sqrt(2) * cos(th) + 2;
y2 = 2*sqrt(2) * sin(th) + 2;
Y(1,:) = x2;
Y(2,:) = y2;

plot(x1, y1);
hold on
plot(x2, y2);
hold off 

for i = 1:r
    for j = 1:c
        A(i,j) = log(norm(X(:,i) - Y(:,j)));
    end
end
        

