

function [v] = prolate_ev(prolate_dat, prolate_ids, xx)
%
% Evaluates the prolate functions.
%
% Input:
%   * prolate_dat : precomputed data structure (prolate_crea).
%   * prolate_ids : which prolates to compute. vector of ids between 0 and prolate_dat.num_prols-1.
%   * xx : vector of points in the interval [0,1] where the prolates should be evaluated
% Output:
%   * v : matrix of evaluate prolates. 
%       each column refers to a different prolate, each row to a different coordinate.
%
%

    assert( prolate_dat.type == 2 )
    
    v = prolate_ZernikeNorm_ex(prolate_dat.p,prolate_dat.N, prolate_dat.cfs(:,prolate_ids+1), xx) ;

end

