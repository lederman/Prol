
function dvec = prolate_xdZernikeNorm_coef(p,N,vec)
%
% Computes the expansion of xf'(x) in the basis of Zernike polynomials,
% where f(x) is given in the basis of Zernike polynomials.
%
% Input:
%  * p,N : Prolate/Zernike parameters.
%  * vec : vector (or multiple vectors in multiple columns) of the
%  coefficients of funtions, expanded in the basis of Zernike polynomials.
%
% Output:
%  * dvec : vector (or multiple vectors in multiple columns) of the
%  coefficients of the expansion of xf'(x). 
%


    dvec = 0 * vec;
    tmpvec = vec;
    tmpjacvec = 0*vec;
    
    for n=size(vec,1)-1:-1:0
        dvec(n+1,:) = dvec(n+1,:) + sqrt(2*n + N + p/2 + 1)/sqrt(2)/(n + N + p/2 + 1) * tmpjacvec(n+1,:);
        dvec(n+1,:) = dvec(n+1,:) + (2*n+N) * (tmpvec(n+1,:));
        
        if (n==0)
            break
        end
        tmpjacvec(n,:)= (n + N + p/2)/(n + N + p/2 + 1) * tmpjacvec(n+1,:);
        tmpjacvec(n,:) = tmpjacvec(n,:) + (2*(n + N) + p)* sqrt(2*(2*n + N + p/2 + 1)) * tmpvec(n+1,:); % Jacobi component        
    end
    

end

