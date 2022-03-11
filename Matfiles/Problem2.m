clc;
clear all;

N=[10, 20, 50, 100, 200, 500, 1000];%, 2000];%, 5000];
%N=[5];
%Matrix Construction
for j=1:length(N)
n=N(j);
x=(0:n-1)./(n-1);
y=2+(0:n-1)./(n-1);
A=log(abs(repmat(x',1,n)-repmat(y,n,1)));
for k=1:10
%Calculation of Lagrange Interpolation
cx=0.5+0.5*cos((2*(1:k)-1)*pi/(2*k));%Interval [0 1]
cx=sort(cx);
cy=2.5+0.5*cos((2*(1:k)-1)*pi/(2*k));%interval [2 3]
cy=sort(cy);
Kxy=log(abs(repmat(cx',1,k)-repmat(cy,k,1)));
L=[];
for i=1:k
L=[L lagranL(x,cx,i)];
end
U=[];
for i=1:k
U=[U lagranL(y,cy,i)];
end
At=L*Kxy*U';
err(j,k)=norm(A-At,2)/norm(A);
end
r(j)=rank(A);
end

function L=lagranL(x,c,i)
n=length(x);
r=length(c)-1;
if r==0
    L=ones(n,1);
else
C=c(i);
c(i)=[];
t1=repmat(x',1,r)-repmat(c,n,1);
t2=prod(C-c);
L=prod(t1,2)./t2;
end
end