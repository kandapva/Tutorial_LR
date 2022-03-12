function [U,V]=partial_ACA(A,mu)
%   This routine computes low rank decomposition of A \approx UV
%   A -  size mxn 
%   U,V' - size mxr 
%   mu - default to pow(10,-10)
if ~exist('mu','var')
    mu=1e-10;
end
k=1;
[m,n]=size(A);
% i,j holds the index of row and column that has max pivot
i=zeros(m,1);
j=zeros(n,1);
i(k)=1;Z=[];
vq=A(i(k),:);
% Initialize U,Vrank 1 approximation 
while norm(vq)==0 && k<m
    Z(end+1)=i(k);
    k=k+1;
    i(k)=k;
    vq=A(i(k),:);
end
if k<m
    j(k)=find_max(vq);
    piv=1/vq(j(k));
    V(k,:)=piv*vq;
    U(:,k)=A(:,j(k));
end
k=k+1;
err=norm(U(:,k-1))*norm(V(k-1,:));
% The general routine that constructs rank 1 approximation at each
% iteration
while err>mu && k<m
   vq=zeros(size(vq));
   i(k)=find_max(U(:,k-1),Z);
   for l=1:k-1
        vq=vq-U(i(k),l)*V(l,:);
   end
   vq=vq+A(i(k),:);
   % The err computed at each iteration is O(n) in opposed to O(n^2) in
   % Cross Approximation
 if norm(vq)~=0
    Z(end+1)=i(k);
    j(k)=find_max(vq);
    piv=1/vq(j(k));
    V(k,:)=piv*vq;
    U(:,k)=A(:,j(k));
    for l=1:k-1
        U(:,k)=U(:,k)-U(:,l)*V(l,j(k));
    end
   err = norm(U(:,k-1))*norm(V(k-1,:));%err=norm(A-U*V);
   k=k+1;
 else
     Z(end+1)=i(k);
     err = norm(U(:,k-1))*norm(V(k-1,:));%err=norm(A-U*V);
 end
end
end
% Finds the index of the maximum element in the given vector
function [maxi]=find_max(u,z)
u=abs(u);
ind=find(u==max(u));
if ~exist('z','var')
    z=[];
end
if isempty(z)
    maxi=ind(1);
else
    for i=1:length(ind)
       if  isempty(find(z==ind(i)))
           maxi=ind(i);break;
       else
           maxi=0;
       end
    end
end
if maxi==0
    u(ind)=0;
    maxi=find_max(u,z);
end  
end