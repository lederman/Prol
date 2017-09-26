function v = prolate_JacobiP_ex(a,b,cfsvec,xx)   
%
% Evaluates functions expanded in Jacobi Polynomials  P^{a,b}_n(x) (not normalized)
%
%  v(i,j) = \sum_{q=0}^{k-1} cfsvec(q+1,j) P^{a,b}_q(x_i)
%
% Input:
%   * a,b : the a,b parameters of the Jacobi polynomials to use here.
%   * cfsvec : k x m matrix. 
%       Columns of coefficients, each column has the k coefficients of an expansion 
%       in Jacobi polynoimals of order k-1 for one of the m different functions.
%   * xx :  a vector of length l. 
%       each entry is a value of x where each one of the k expansions should be evaluated.
% Output:
%   * v : l x m matrix.
%       The j-th column is the j-th function evaluated at the l points.
%
%

    xx = xx(:); 

    jcm1 = 0*xx(:)+1;    
    jc0  = 1/2*(a-b+(2+a+b)*xx(:));
    jcp1 = 0*xx(:);    
    if (size(cfsvec,1)>0)
        v= jcm1 *cfsvec(1,:); 
    else
        v=zeros([length(xx),size(cfsvec,2)],class(xx));
    end
    if size(cfsvec,1)>1
        v=v+jc0 * cfsvec(2,:);
        for k=1:size(cfsvec,1)-2
            jcp1(:) = (((2*k+a+b+1)*(a^2-b^2)+(2*k+a+b)*(2*k+a+b+1)*(2*k+a+b+2)*xx(:)).*jc0 -2*(k+a)*(k+b)*(2*k+a+b+2)*jcm1) / (2*(k+1)*(k+a+b+1)*(2*k+a+b));
            jcm1(:) = jc0(:);
            jc0(:)  = jcp1(:);
            v=v+jc0*cfsvec(k+2,:);
        end
    end
end
