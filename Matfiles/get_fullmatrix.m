function A = get_fullmatrix(X)
% Naive Implementation 
N =length(X);
A = zeros(N,N);
for i=1:N
    for j=1:N
        A(i,j) = X(i) - X(j);
    end
end

%A = ri-rj;
end