function [U,V]=cross_app(A,tol)
error =1;
iter=1;
M=A;
U=[];
V=[];
if ~exist('tol','var')
    tol=1e-12;
end
while error>tol 
[is,js]=maxmat(abs(A));
del=A(is,js)
a=A(:,js)/del;
U(:,iter)=a;
b=A(is,:);
V(iter,:)=b;
M=M-a*b;
A=A-a*b
error=norm(M);
iter=iter+1;
if iter>8
    break;
end
end
end
function [ a,b ] = maxmat( c )
as=size(c);
total_ele=numel(c);
[~,I]=max(c(:));
r=rem(I,as(1));
a=r;
b=((I-a)/as(1))+1;
if a==0
    a=as(1);
    b=b-1;
else
    a=r;
    b=b;
end
end