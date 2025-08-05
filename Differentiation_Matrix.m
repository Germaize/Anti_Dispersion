function [ D,x ] = Differentiation_Matrix(N)

%   
%   [ 2 -1  1 -1  1 -1 ... 2  ]'
%
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';

% chebyshev points -> vector
n = (0:N)'     ;
x = cos(pi*n/N); 

% replicate vector N+1 times
X = repmat(x,1,N+1);

% in  column j : we have x_i - x_j
dX = X-X';

%  c*(1./c)' =   c_i/c_j
%  dX+(eye(N+1)) = x_i - x_j +1
D  = (c*(1./c)')./(dX+(eye(N+1)));  

D  = D - diag(sum(D'));   

% More details in Quarteroni & Zang 

end

