function [ratios , ratios_consecutive] = prolate_crea_eigRatios(prolate_dat)
%
% Compute the ratios between consecutive eigenvlues: 
%    \gamma_{n+1} / \gamma_{n}, 
% and the ratios between these eigenvalues and the first eigenvalue as
% reference:
%   \gamma_{n} / \gamma_{0}
%
% Input:
%  * prolate_dat : the prolates datastructure created by prolate_crea.
% Ourput:
%  * ratios: Ratios between the eigenvalues and the 0-th eignevalue:
%       \gamma_{n} / \gamma_{0}
%       note that this vector is zero-bases, so that:
%       ratios(1) = 1 , ratios(n) = \gamma_{n+1} / \gamma_{0}
%  * ratios_consecutive : Ratios between consecutive eigenvalues.
%       again, zero based, so the first ration is meaningless. 
%

    assert( prolate_dat.type == 2 ) 

    % get the x*d/dx of the prolates, expanded in the same Zernike
    % polynomials expansion.
    dvec = prolate_xdZernikeNorm_coef(prolate_dat.p, prolate_dat.N, prolate_dat.cfs);
    
    % Get the numerator and denominator of of relations between eigenvalues
    % through inner products between the espansions.
    int_num = sum(circshift(dvec,[0,1]).*prolate_dat.cfs,1); % int_p1(k) = \int_0^1 x \Psi_{k-1}'(x) \Psi_k(x) x^(p+1) dx
    int_den = sum(dvec.*circshift(prolate_dat.cfs,[0,1]),1); % int_m1(k) = \int_0^1 x \Psi_{k}'(x) \Psi_{k+1}(x) x^(p+1) dx
    
    % Compute the consecutive ratios:
    ratios_consecutive = int_num./(int_den);   %
    ratios_consecutive(1) = 1;
    
    % Cumulative ratio
    ratios = cumprod(ratios_consecutive);

end
