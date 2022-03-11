% This Code Tests the partially pivoted ACA
clc;
clear all;
N=[10,50,100,500,1000,2000];
k=1;
err = zeros(length(N),7);
a_time = zeros(length(N),7);
for i=1:length(N)
    for j=1:7
        A=kernel_matrix(N(i),j,8);
        tic;
        [u,v]=partial_ACA(A,1e-10);
        a_time(i,j) =toc();
        err(i,j)=norm(A-u*v,2)/norm(A,2);
    end
end