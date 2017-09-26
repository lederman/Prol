
function [vdiag, voffdiag] = prolate_diffop_mat_tridiag(c,p,N,maxk)  
%
% Computes the matrix representation of the differential operator,
% in the basis of Zernike polynomials. 
%
% Input:
%   * c,p,N : prolate parameters.
%   * maxk  : matrix truncations: the dimensionality of the matrix is k+1
% Output:
%   * vdiag : vector of the elements of the matrix diagonal
%   * voffdiad : vector of coefficients of the matrix off diagonal
%
% Note that the matrix is symmetric tridiagonal. All other elements are
% zeros.
%


    vdiag = zeros(maxk+1,1);
    voffdiag = zeros(maxk,1);
    
    for k=0:maxk
        vdiag(k+1,1) = prolate_diffop_mat_diag_element(c,p,N,k);
    end
    for k=0:maxk-1
        voffdiag(k+1,1) = prolate_diffop_mat_offdiag_element(c,p,N,k);
    end
    
end




function v = prolate_diffop_mat_offdiag_element(c,p,N,k)
%
% helper function for prolate_diffop_mat_tridiag.
% offdiagonal elements.
%
    nn = k+1;
    NN = N+p/2;
    if nn<=0 
        v=0;
    else
        v = -(c^2*nn)/((2*nn+NN)*(2*nn+NN+1)) * (nn+NN)/(sqrt(1-2/(1+2*nn+NN)));
    end    
end



function v = prolate_diffop_mat_diag_element(c,p,N,k)
%
% helper function for prolate_diffop_mat_tridiag.
% elements on the diagonal.
%
    NN=N+p/2;
    if (NN==0)&&(k==0)
        v = -( prolate_diffop_mat_kappa(p,N,k) +c^2/2 );
    else
        v = -( prolate_diffop_mat_kappa(p,N,k) +c^2*(2*k*(k+1)+NN*(2*k+NN+1))/((2*k+NN)*(2*k+NN+2)) );
    end
end




function kap = prolate_diffop_mat_kappa(p,N,k)
%
% helper function for prolate_diffop_mat_tridiag.
%
    NN = N+p/2;
    kap = (NN+2*k+1/2)*(NN+2*k+3/2);
end


