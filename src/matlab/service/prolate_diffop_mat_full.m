
function [mat, vdiag, voffdiag ] = prolate_diffop_mat_full(c,p,N,maxk)
%
% Computes the full matrix representation of the differential operator,
% in the basis of Zernike polynomials. 
%
% Input:
%   * c,p,N : prolate parameters.
%   * maxk  : matrix truncations: the dimensionality of the matrix is k+1
% Output:
%   * mat : the matrix
%   * vdiag : vector of the elements of the matrix diagonal
%   * voffdiad : vector of coefficients of the matrix off diagonal
%
% Note that the matrix is symmetric tridiagonal. All other elements are
% zeros.
%

    mat = zeros(maxk+1);
    [vdiag, voffdiag] = prolate_diffop_mat_tridiag(c,p,N,maxk);
    
    for k=0:maxk
        mat(k+1,k+1) = vdiag(k+1);        
    end
    for k=0:maxk-1
        mat(k+1,k+2) = voffdiag(k+1);        
        mat(k+2,k+1) = voffdiag(k+1);
    end
    
end

