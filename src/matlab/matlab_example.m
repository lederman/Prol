function matlab_example()
%
% Example for using the matlab version of this code. 
%

    % add to path
    matlab_addpath_prol_src();
    
    % prolate parameters
    c=pi*20;   % band limit
    D=3;       % prolate on 3-D ball    
    N = 1;     % prolate angular frequency parameter: N=0,1,.... , see notes. 

    minEigenvalRatio = 10^-40; % the truncation parameter. magnitude of the smallest eigenvalue that we would like to keep
    
    % technical parameters:
    matdim = 1000;                        % size of the matrix representation of the differential operator, see notes. If this is too small, a warning should appear. This will be removed in future versions.    
    prolate_crea_options.isfixfirst = 1; % this option should be set to 1 for more accurate eigenvalues.
    
    %
    %
    %
    disp('Precomputation')
    tic
    [prolate_dat, iserr , ~] = prolate_crea(c,D,N,minEigenvalRatio, matdim, prolate_crea_options);    
    toc
    if (iserr~=0)
        disp('An error has occured!');
    end
    
    h0 = figure;
    semilogy([0:prolate_dat.num_prols-1],(abs(prolate_dat.nu)),'LineWidth',3)   
    ylim([10^-30,3])
    xlabel('n')    
    set(gca,'FontSize', 12);
    ylabel('|\nu_n|','FontSize', 14)
    
    %
    % Additional eigenvalues (see documentation):
    %
    % The other versions of the integral operator eigenvalues are available
    % in prolate_dat.gam , prolate_dat.bet, prolate_dat.alp.
    %
    % The eigenvalues of the differential operator are availavle in 
    % prolate_dat.chi
    %
    
    %
    %
    %
    disp('Evaulate radial part at given points')
    
    prol_ids = [0:prolate_dat.num_prols-1];  % ids of the prolates we wish to evalualte. Note that the index of the first prolate is 0. We evaluate all of the, and then display a selection below.
    xx = linspace(0,1,1000);                 % points where the prolate should be evaluated. Must be in [0,1].
    
    tic
    [v] = prolate_ev(prolate_dat, prol_ids , xx);
    toc   
    
    h1 = figure; 
    funcid = [0:1,2,5,10];
    funcid( funcid > prolate_dat.num_prols-1 ) = []; % remove indices that are out of range, in case the user changes the parameters of this example.
    plot(xx,v(:,funcid+1),'LineWidth',2);  % +1 because the prolate indices start at 0
    xlabel('x')
    lgd=legend(num2str(funcid'));
    set(gca,'FontSize', 12);
    ylabel('\Phi_{N,n}(x)','FontSize', 14)    
    


end