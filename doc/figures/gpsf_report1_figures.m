%
% prol (Matlab(R) version)
% Demosntration code for computing generalized prolate spheroidal functions.
% 
% Author: Roy R. Lederman
% http://roy.lederman.name
% http://github.com/lederman/prol
%
% This code generates the figures for the paper gpsf_report1.tex 
%



function prol_report_figures()

    % run matlab_addpath_prol_src() in /src/matlab before running this code. 


    report_part001()

    
   
end



function report_part001()

    %
    % eigenvalues
    %
    c=pi*20;
    D=3;
    h1=figure;
    h2=figure;
    matdim = 800;
    minEigenvalRatio = 10^-40;
    prolate_crea_options.isfixfirst = 1; 
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
    %
    %
    h1=figure;
    h2=figure;    
    prolate_crea_options.isfixfirst = 1; 
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
    %
    %
    c= 10 * pi;
    D=3;
    N=1;
    matdim = 800;
    xx = linspace(0,1,1000);
    minEigenvalRatio = 10^-30;
    
     
    
    tic
    prolate_crea_options.isfast = 0; 
    [prolate_dat, iserr , ~] = prolate_crea(c,D,N,minEigenvalRatio, matdim, prolate_crea_options);
    toc
    
    tic
    [v] = prolate_ev(prolate_dat, [0:prolate_dat.num_prols-1], xx);
    toc
    
    figure; plot(xx,v(:,1));
    figure; plot(xx,v(:,[1:3:end]));
    
    
end



