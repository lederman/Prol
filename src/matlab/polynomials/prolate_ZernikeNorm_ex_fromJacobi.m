
function v = prolate_ZernikeNorm_ex_fromJacobi(p,N,cfsvec,xx) 
%
% Evaluates functions expanded in the basis of normalized Zernike polynomials
% using Jacobi polynomials.
%
%  v(i,j) = \sum_{q=0}^{k-1} cfsvec(q,j) \hat{R}_{N,n,p}_q(x_i)
%
%  Using Jacobi polynomials:
%   \hat{R}_{N,n,p}_q(x_i) = (-1)^n \sqrt{2(2n+N+p/2+1)} x^N P^{(N+p/2,0)}_n(1-2x^2)
%
% Input:
%   * p,N : the p,N parameters of the Zernike polynomials to use here.
%   * cfsvec : k x m matrix. 
%       Columns of coefficients, each column has the k coefficients of an expansion 
%       in Zernike polynoimals of order k-1 for one of the m different functions.
%   * xx :  a vector of length l. 
%       each entry is a value of x where each one of the k expansions should be evaluated.
% Output:
%   * v : l x m matrix.
%       The j-th column is the j-th function evaluated at the l points.
%
% 
    
    b=0;
    a=N+p/2;
    yy = 1-2 * xx(:).^2;
    
    K = size(cfsvec,1)-1;
    
    cfsvec_jac = cfsvec;
    cfsvec_jac = bsxfun(@times, (-1).^[0:K]', cfsvec_jac);
    cfsvec_jac = bsxfun(@times, sqrt(2*(2*[0:K]' + N + p/2 + 1)), cfsvec_jac);
    
    v = prolate_JacobiP_ex(a,b,cfsvec_jac,yy) ;
    
    v=bsxfun(@times, xx(:).^N , v);
 
end
