function [U,S,V] = rand_svd(A,k)
[~,n] = size(A);
%%% rand_svd takes the matrix A and decomposes it into U,S,V
%%% Given matrix A with size m x n, U is m x k,V is n x k, S is kxk
%%% we use a constant over-sampling parameter p = 5 such that r = k+p
p = 5;
r = k + p;
%%% STEP 1: Range Finder, find a Q(m x r) such that ||A - QQ'A|| < \epsilon
Rmat = randn(n,r);
Q = A*Rmat; % Size of Q is m x rCost is mnr
[Q, ~] = qr(Q);% Size of Q is m x r m*r^2
%%% STEP 2: 
A = Q'*A; % Size of A is rxn
[u,S,V] = svd(A,"econ");%nr^2
U         = Q*u(:,1:k);
S         = S(1:k,1:k);
V         = V(:,1:k);
end