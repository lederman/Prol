
function v = prolate_ZernikeNorm_ex(p,N,cfsvec,xx) 
%
%
% Evaluates functions expanded in the basis of normalized Zernike
% polynomials.
%
%  v(i,j) = \sum_{q=0}^{k-1} cfsvec(q,j) \hat{R}_{N,n,p}_q(x_i)
%
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

    v = prolate_ZernikeNorm_ex_fromJacobi(p,N,cfsvec,xx) ;
    
    
end
