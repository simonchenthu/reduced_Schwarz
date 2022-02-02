function beta = legendre_recurrence(N)
% recurrence -- recurrence coefficients for Legendre polynomials
%
%     Calculates coefficients of the first N+1 recurrence formulas 
%     for the (normalized) Legendre polynomials.
%
%     (n+1)/sqrt((2n+1)*(2n+3))*P_{n+1}+n/sqrt((2n+1)*(2n-1))*P_{n-1}=xP_n,
%     n=0,1,...,N


beta = [0:1:N];
beta = (beta+1)./sqrt((2*beta+1).*(2*beta+3));

return