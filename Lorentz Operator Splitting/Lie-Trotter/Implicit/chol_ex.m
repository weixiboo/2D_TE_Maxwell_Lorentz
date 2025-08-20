
N = 10000;

A = diag(ones(N,1)*2) + diag(ones(N-1,1)*-1,-1) + diag(ones(N-1,1)*-1,1);
b = ones(N,1);
R = chol(A);

tic;
x = A\b;
toc

tic;
opts_1.LT = true;
x_chol = linsolve(R',b,opts_1);

opts_2.UT = true;
x_chol = linsolve(R,x_chol,opts_2);
toc

max(abs(x-x_chol))
