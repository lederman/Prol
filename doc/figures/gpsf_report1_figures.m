%
% prol 
% Demosntration code for computing generalized prolate spheroidal functions.
% (Matlab(R) version)
% 
% Author: Roy R. Lederman
% http://roy.lederman.name/
% http://github.com/lederman/prol
%
% This code generates the figures for the paper gpsf_report1.tex 
%



function gpsf_report1_figures()

    % run matlab_addpath_prol_src() in /src/matlab before running this code. 

    file_header = 'gpsf_report1_'
    report_part001(file_header)

    
   
end



function report_part001(file_header)

    %
    % sample eigenvalues figures
    %
    c=pi*20;
    D=3;
    h1=figure; % eigenvalues magnitude |\nu|
    h2=figure; % eigenvalues magnitude is close to one: |1-|\nu||
    matdim = 800;
    minEigenvalRatio = 10^-40;
    prolate_crea_options.isfixfirst = 1; 
    Ns = [0:5:20];
    for j1=1:length(Ns)
        N=Ns(j1);
        tic
        [prolate_dat, iserr , ~] = prolate_crea(c,D,N,minEigenvalRatio, matdim, prolate_crea_options);
        toc
        figure(h1)        
        semilogy([0:prolate_dat.num_prols-1],(abs(prolate_dat.nu)),'LineWidth',3)
        hold on
        figure(h2)
        semilogy([0:prolate_dat.num_prols-1],max(abs(1-abs(prolate_dat.nu)),10^-20),'LineWidth',3)
        % the max is taken to avoid log(0) when |\nu| = 1 exactly.
        hold on
    end
    figure(h1)
    ylim([10^-30,3])
    xlabel('n')
    lgd=legend(num2str(Ns'));
    %title(lgd,'N=')
    set(gca,'FontSize', 12);
    ylabel('|\nu_n|','FontSize', 14)
    
    print([file_header,'D3_eigenvals.png'],'-dpng')
    
    
    figure(h2)
    ylim([10^-16,3])
    xlabel('n')
    lgd=legend(num2str(Ns'));
    %title(lgd,'N=')
    set(gca,'FontSize', 12);
    ylabel('|1-|\nu_n||','FontSize', 14)
    
    print([file_header,'D3_eigenvals_to_one.png'],'-dpng')

    
    
    %
    % Without fixing the coefficients of the first eigenvector
    %
    h1=figure;
    h2=figure;    
    prolate_crea_options.isfixfirst = 0; 
    Ns = [0:3:20];
    for j1=1:length(Ns)
        N=Ns(j1);
        tic
        [prolate_dat, iserr , ~] = prolate_crea(c,D,N,minEigenvalRatio, matdim, prolate_crea_options);
        toc
        figure(h1)
        semilogy([0:prolate_dat.num_prols-1],(abs(prolate_dat.nu)),'LineWidth',3)
        hold on
        figure(h2)
        semilogy([0:prolate_dat.num_prols-1],abs(1-abs(prolate_dat.nu)),'LineWidth',3)
        hold on
    end
    figure(h1)
    ylim([10^-30,3])
    xlabel('n')
    lgd=legend(num2str(Ns'));
    %title(lgd,'N=')
    set(gca,'FontSize', 12);
    ylabel('\nu_n','FontSize', 14)
    
    
    
    %
    % sample eigenfunctions figures
    %
    c= 20 * pi;
    D=3;
    N=0;
    matdim = 800;
    xx = linspace(0,1,1000);
    minEigenvalRatio = 10^-30;
    prolate_crea_options.isfixfirst = 1; 
    
    
    tic
    [prolate_dat, iserr , ~] = prolate_crea(c,D,N,minEigenvalRatio, matdim, prolate_crea_options);
    toc   
    tic
    [v] = prolate_ev(prolate_dat, [0:prolate_dat.num_prols-1], xx);
    toc    
    h1 = figure; 
    funcid = [0:1,2,5,10];
    plot(xx,v(:,funcid+1),'LineWidth',2);
    xlabel('x')
    lgd=legend(num2str(funcid'));
    set(gca,'FontSize', 12);
    ylabel('\Phi_{N,n}(x)','FontSize', 14)    
    
    print([file_header,'D3_N0_eigenfuncs.png'],'-dpng')
    
    %
    %
    %
    N=1;
    tic
    [prolate_dat, iserr , ~] = prolate_crea(c,D,N,minEigenvalRatio, matdim, prolate_crea_options);
    toc   
    tic
    [v] = prolate_ev(prolate_dat, [0:prolate_dat.num_prols-1], xx);
    toc    
    h1 = figure; 
    funcid = [0:1,2,5,10];
    plot(xx,v(:,funcid+1),'LineWidth',2);
    xlabel('x')
    lgd=legend(num2str(funcid'));
    set(gca,'FontSize', 12);
    ylabel('\Phi_{N,n}(x)','FontSize', 14)    
    
    print([file_header,'D3_N1_eigenfuncs.png'],'-dpng')
    
    
end



