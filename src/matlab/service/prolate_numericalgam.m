
function gam = prolate_numericalgam(prolate_dat, n)
%
% Numerical computation of the n-th eigenvalue of the integral operator.
% Uses the data structure prolate_dat created by prolate_crea.
%
% Generally speaking, this function should only be used for computing the 
% eigenvalue for n=0 by prolate_crea, and should not be used otherwise.
%
% Input:
%  * prolate_dat : precomputed prolate information.
%  * n : the id of the eigenvalue to be computed. 
%       This would usually be n=0.
% Output:
%  * gam : the eigenvalue \gamma_n
%

% Todo: remove matlab dependency

    assert( prolate_dat.type == 2 )

    % find a large enough point
    % TODO: replace with approximate maximum using WKB and/or Newton method search.
    xx0 = linspace(0,1,1000)';
    % prolate
    [v] = prolate_ev(prolate_dat, [n], xx0);
    % weighted prolate: \phi_n(x) = x^{(p+1)/2} \Phi_n (x);
    v=bsxfun(@times, xx0.^((prolate_dat.p+1)/2) , v);
    % find max
    [xmax_v,xmax_id] = max(abs(v));
    xmax_v = v(xmax_id);
    xmax = xx0(xmax_id);
    
    % truncate the vector of coefficients
    vec = prolate_dat.cfs(:,n+1);
    tmpkeep = find(abs(vec) >= prolate_dat.evparam.cfs_eps);
    idskeep=tmpkeep(end);
    vec((idskeep+1):end) = [];
    
    %
    % numerical integration:
    %    
    
    % function to integrate
    fun = @(y) reshape( besselj(prolate_dat.N+prolate_dat.p/2, prolate_dat.c*xmax*y(:)).*sqrt(prolate_dat.c *xmax *y(:)) .*y(:).^((prolate_dat.p+1)/2) .* prolate_ZernikeNorm_ex(prolate_dat.p,prolate_dat.N, vec, y(:)) , size(y) );
    % integration:
    %q1 = integral( fun,0,1 ); % matlab only
    q1 = quad( fun,0,1, eps(xmax_v)*2 ); % compatible with Octave
    
    %
    % The eigenvalue is the ratio:
    %
    gam = q1 / xmax_v;
    

end

