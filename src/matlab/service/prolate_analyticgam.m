
function gam = prolate_analyticgam(prolate_dat, n)
%
% Computation of the n-th eigenvalue of the integral operator.
% Uses the data structure prolate_dat created by prolate_crea.
%
% Generally speaking, this function should only be used for computing the 
% eigenvalue for n=0 by prolate_crea.
%
% Input:
%  * prolate_dat : precomputed prolate information.
%  * n : the id of the eigenvalue to be computed. 
%       This would usually be n=0.
% Output:
%  * gam : the eigenvalue \gamma_n
%

% Todo: remove dependency on undocumented properties of the eigenvector
% computation in matlab.

%
% Note: this function can be more sensitive to the truncation of the list
% of coefficients.
%

    % coefficients of the chosen prolate
    cfs = prolate_dat.cfs(:,n+1);
    
    % extract parameters
    N=prolate_dat.N;
    p=prolate_dat.p;
    c=prolate_dat.c;
    k=[0:length(cfs)-1]';
    
    % parts of the computation
    cfs1  = (-1).^k .* sqrt(2+4*k+2*N+p) .* (2+2*N+p) /2;
    cfs2n = (N+p/2 + k);
    cfs2n(1) = gamma(1+ N+p/2 );
    cfs2d = k;
    cfs2d(1) = 1;
    cfs2 = cumprod(cfs2n./cfs2d);
    
    %
    % Safety truncation to avoid inf.
    %
    mytrunc1 = find(abs(cfs2)>realmax*10^-10);
    mytrunc2 = find(abs(cfs)<realmin*10^10);
    mytrunc = min([mytrunc1,mytrunc2]);
    if ~isempty(mytrunc)
        cfs = cfs(1:mytrunc);
        cfs1 = cfs1(1:mytrunc);
        cfs2 = cfs2(1:mytrunc);        
    end
    
    %
    % Compute \gamma
    %
    num = 2^(-(N+p/2+1)) * c^(N+p/2+0.5) * sqrt(2+2*N+p) * cfs(1);
    denom = sum( cfs.*cfs1.*cfs2 );
    
    gam = num/denom;
  
end

