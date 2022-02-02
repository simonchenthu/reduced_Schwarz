function [X,W] = legendre_quad(n) 

% Golub-Welsch algorithm for Gauss-Legendre quadrature with n+1 nodes in [-1,1]
% X: zero points of (n+1)th order Legendre polynomials
% W: Gauss-Legendre weights corresponding to zero points

beta = legendre_recurrence(n);

A = diag(beta(1:end-1), 1) + diag(beta(1:end-1), -1);

[V,X]=eig(A);
X = diag(X);
W = 2*V(1,:)'.^2;    % w_j=v(1,j)^2 * int_{-1}^{1}rho(x)dx 


return