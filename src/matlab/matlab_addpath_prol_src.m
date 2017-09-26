%
% Add to path
%
function matlab_addpath_prol_src()
    path_to_pkg = fileparts(mfilename('fullpath'));

    addp = @(d)(addpath(fullfile(path_to_pkg, d)));

    addp('');
    addp('polynomials');
    addp('service');

end
