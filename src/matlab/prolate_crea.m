
function [prolate_dat, iserr , prolate_dat_tmp] = prolate_crea(c, D, N, minEigenvalRatio, matdim , prolate_crea_options)
%
% prolate_crea creates a data-structure for computing a family of
% generalized prolate spheroidal functions for dimension D and order N.
%
% Input:
%  * c : prolate truncation frequency
%  * D : prolate dimension (D=p+2)
%  * N : prolate order
%  * minEigenvalRatio : keep only prolates the with eigenvalue \gamma s.t.
%       | \gamma_n | >  c^{-1/2} * minEigenvalRatio .
%       The reason is that for small n and large c, c^{1/2}| \gamma_n | ->  1
%  * matdim : the dimensionality of teh matrix used in precomputation. 
%       This is a technical parameter which will be removed in future versions. 
%       If this number is too small, a warning will be generated.
%       If this number is too large, the precomputation can be slow.
%  * prolate_crea_options : optional patameters
%       isfast : run faster computation without fixing the eigenvectors? 
%                This will be removed in future versions.
%
% Output:
%  * prolate_dat : data structure to be used in prolate_ev
%  * iserr    : error code
%               0 : no error
%               1 : empty set of prolates.
%               10 : matdim may be too small (based on number of prolates kept)
%               100 : matdim may be too small (based on number of coefficients kept)
%
%

%
% TODO: 
%  * Remove Matlab eig for more accurate computation of eigenvectors'
%    elements without the second pass on the eigenvectors. 
%  * Introduce wrapper to compute matdim
%  * Introduce warpper to compute for all N. 
%  * Add user control of accuracy.
%  

    assert(round(D)==D)     % Integer dimension
    assert(D>1)             % One dimensional case to be treated separately 
    assert(c>0)             % Physical c
    assert(N>=0)            % Physical N
    assert(round(N)==N)
    
    assert(minEigenvalRatio<1)      % otherwise, what is the point?    
    assert(minEigenvalRatio > 0)    % Eigenvalues must be truncated
    assert(matdim > 10)             % The matrix for computing the coefficients cannot be too small
    
    isfast = 1; % don't bypass the eigenvector correction
    isfixfirst = 1; % don't bypass the eigenvector correction
    if exist('prolate_crea_options')
        if isfield(prolate_crea_options,'isfast')
            isfast = prolate_crea_options.isfast;
        end
        if isfield(prolate_crea_options,'isfixfirst')
            isfixfirst = prolate_crea_options.isfixfirst;
        end
    end
        
    
    
    %
    % Parameters
    %
    iserr = 0;
    prolate_dat.type = 2; 
    prolate_dat.c = c;
    prolate_dat.D = D;
    prolate_dat.p = D-2;
    prolate_dat.N = N;
    prolate_dat.creaparam.minEigenvalRatio = minEigenvalRatio;
    prolate_dat.creaparam.matdim = matdim;
        
    prolate_dat.evparam.cfs_eps = eps(1.0)/100;
    
    %
    % differential equation eigenproblem in matrix form, the eigenvectors
    % ate the coefficients of the prolats.
    %    
    
    % TODO: replace the full matrix operation.
    [mat, vdiag, voffdiag] = prolate_diffop_mat_full(prolate_dat.c, prolate_dat.p ,prolate_dat.N , prolate_dat.creaparam.matdim-1);
        
    [u,d] = eig(mat);
    % Note that Matlab eig truncates some small coefficients by setting them to 0. 
    [eigvals,eigvals_order] = sort(diag(d),'descend');
    eigvecs = u(:,eigvals_order);
    
    
    eigvecs = bsxfun(@times, eigvecs, sign(eigvecs(1,:)).*(-1).^[0:matdim-1] ); % standard sign. Assumes accurate first element.
    
    %
    % temporary data structure that stores data before truncation
    %
    prolate_dat_tmp = prolate_dat;
    prolate_dat_tmp.cfs = eigvecs;
    prolate_dat_tmp.diffeigs = eigvals;
    
    
    % fix the first eigenvector to reduce the scope of numerical inaccuracy
    % due to eigenvector truncation in Matlab's eig.
    if (isfixfirst==1)
        prolate_dat_tmp.cfs(:,1) = prolate_crea_fix_eigenvec(mat, prolate_dat_tmp , 1);
    end
    
    % compute the first eigenvalue        
    %gam0num = prolate_numericalgam(prolate_dat_tmp, 0);% TODO: replace with approximate maximum using WKB and/or Newton method search.    
    gam0 = prolate_analyticgam(prolate_dat_tmp, 0);         
    % Compute the rest of the eigenvalues through recurrsion    
    [ratios , ~] = prolate_crea_eigRatios(prolate_dat_tmp);    
    prolate_dat_tmp.gams = gam0 * ratios;
    prolate_dat_tmp.nu_abs = abs(gam0 * ratios * prolate_dat.c^(1/2));        
    % Find where to truncate
    ids_prolate_to_discard = find( prolate_dat_tmp.nu_abs <= prolate_dat.creaparam.minEigenvalRatio );
    ids_prolate_to_keep = [1:min(ids_prolate_to_discard)-1+1];
    if (isempty(ids_prolate_to_keep))
        warning('No prolates to keep');
        iserr = 1;
        return
    end
    
    
    % Note that Matlab eig truncates the small coefficients by setting them to 0. 
    abs_cfs = max( abs(prolate_dat_tmp.cfs (:,[1:ids_prolate_to_keep(end)])), [], 2) ;
    cfs_to_keep = find( abs_cfs >= prolate_dat.evparam.cfs_eps );
    cfs_to_keep = [1:max(cfs_to_keep)];
    
    % truncated coefficients
    prolate_dat.cfs = prolate_dat_tmp.cfs( 1:cfs_to_keep(end) , 1:ids_prolate_to_keep(end) ) ;
    
    % the various forms of the integral operator eigenvalues
    prolate_dat.gam = gam0 * ratios(1:ids_prolate_to_keep(end) );
    prolate_dat.bet = prolate_dat.gam/(prolate_dat.c^((prolate_dat.p+1)/2));
    prolate_dat.alp = prolate_dat.bet * (1i)^prolate_dat.N * (2*pi)^(1+prolate_dat.p/2);
    
    prolate_dat.nu = (1i)^prolate_dat.N * prolate_dat.c^(1/2) * prolate_dat.gam;
    
    % the differential operator eigenvalues
    prolate_dat.chi = eigvals(1:ids_prolate_to_keep(end));
    
    prolate_dat.num_prols = ids_prolate_to_keep(end);
    
    %
    % Warnings
    % Take plenty of margin for the truncation.
    %
    if (ids_prolate_to_keep(end)+20 >= matdim)
        warning('prolate_crea: insufficient margin in matrix size (number of prolates)')
        iserr = iserr+10;
    end
    if (cfs_to_keep(end)+40 >= matdim)
        warning('prolate_crea: insufficient margin in matrix size (number of coefficients)')
        iserr = iserr+100;        
    end
    
    
    
end


function v = prolate_crea_fix_eigenvec(mat, prolate_dat , jj)
    
    tmpmat = mat-eye(prolate_dat.creaparam.matdim)*...
        (prolate_dat.diffeigs(jj)+(prolate_dat.diffeigs(jj+1)-prolate_dat.diffeigs(jj))/10^5/(jj+1) );
    v=tmpmat\prolate_dat.cfs(:,jj);
    v = v/norm(v);
    
end