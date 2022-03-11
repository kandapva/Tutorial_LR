function A=kernel_matrix(n,ch,r)
if ~exist('r','var')
    r=1;
end
x=linspace(-1,0,n);
x=repmat(x,[n,1]);
y=linspace(0,1,n)';
y=repmat(y,[1,n]);
d=x-y;
    switch ch
        case 1 
            A=1+(d.^2);
        case 2
            A=sqrt(1+(d.^2));
        case 3
            A=1./(1+(d.^2));
        case 4
            A=1./sqrt(1+(d.^2));
        case 5
            A=exp(-d.^2);
        case 6
            A=exp(-d);
        case 7%Random Matrix
            A=randn(n,r);
            A=A*A';
        case 8%Complex Entries
            A=log(1+d);  
    end  
end