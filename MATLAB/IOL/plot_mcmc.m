function [] = plot_mcmc()

    close all;
    clear all;
    clc;
    set(0,'DefaultFigureVisible','on');

    %%% generate an output directory
    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end

    %%% Vero cells IOL %%%
    data_IOL_PFU_Vero       = readtable('../../Supplementary_Data.xlsx','Sheet','IOL PFU Vero');
    data_IOL_RNA_Vero       = readtable('../../Supplementary_Data.xlsx','Sheet','IOL RNA Vero');
    data_IOL_PFU_Vero_Decay = readtable('../../Supplementary_Data.xlsx','Sheet','IOL PFU Vero Decay');
    data_IOL_RNA_Vero_Decay = readtable('../../Supplementary_Data.xlsx','Sheet','IOL RNA Vero Decay');

    stdev_IOL_PFU_Vero       = std(log10(data_IOL_PFU_Vero(:,2:end)),0,2,'omitnan');
    stdev_IOL_RNA_Vero       = std(log10(data_IOL_RNA_Vero(:,2:end)),0,2,'omitnan');  
    stdev_IOL_PFU_Vero_Decay = std(log10(data_IOL_PFU_Vero_Decay(:,2:end)),0,2,'omitnan');  
    stdev_IOL_RNA_Vero_Decay = std(log10(data_IOL_RNA_Vero_Decay(:,2:end)),0,2,'omitnan');

    stdev = {{stdev_IOL_PFU_Vero_Decay;
            stdev_IOL_RNA_Vero_Decay};
            {stdev_IOL_PFU_Vero;
            stdev_IOL_RNA_Vero}};

    data = {{data_IOL_PFU_Vero_Decay;
            data_IOL_RNA_Vero_Decay;};
            {data_IOL_PFU_Vero;
            data_IOL_RNA_Vero}};

    colors = {[0 0.5 0],[0 0 0];...
              [1 0 1],[0 0 0]};

    YLims = {[-2.6 10.6];...
             [-0.5 10.5]
             };   
    YTicks = {[0,2,4,6,8,10];...
              [0,2,4,6,8,10]};

    figlabs     = {'(A) Viral decay','(B) Infection 0.1 PFU/cell','(C) Practical identifiability analysis'};

    fontsize = 12;

    time_decay  = [0,2,4,8,16,24,48,72,96,120];
    time_infec  = [1.5,3,5,9,17,25,49,73,97,121];
    time = {time_decay;...
            time_infec};

    XLims = {[-4 124];...
             [-3 125]};
    XTick = {[0,8,24,48,72,96,120];...
             [1,9,25,49,73,97,121]};

    T0      = 1;
    MOI     = 0.1;            
    T0_pfu_Vero  = 3e+5;
    V0_pfu_Vero  = MOI*T0_pfu_Vero;

    tmax     = 121;
    tt_decay = 0:0.1:120;
    tt_infec = 0:0.1:121;
    t = {tt_decay;...
         tt_infec};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    outputvec_decay_Vero_pfu     = readmatrix('output_mcmc_IOL_Vero_decay_pfu.txt');
    outputvec_decay_Vero_rna     = readmatrix('output_mcmc_IOL_Vero_decay_rna.txt');
    outputvec_infec_Vero_pfu     = readmatrix('output_mcmc_IOL_Vero_pfu.txt');
    outputvec_infec_Vero_rna     = readmatrix('output_mcmc_IOL_Vero_rna.txt');

    %%% 95% CI %%%
    decay_Vero_pfu_95CI_min = quantile(outputvec_decay_Vero_pfu,0.025,1);
    decay_Vero_pfu_95CI_max = quantile(outputvec_decay_Vero_pfu,0.975,1);
    decay_Vero_rna_95CI_min = quantile(outputvec_decay_Vero_rna,0.025,1);
    decay_Vero_rna_95CI_max = quantile(outputvec_decay_Vero_rna,0.975,1);
    
    infec_Vero_pfu_95CI_min = quantile(outputvec_infec_Vero_pfu,0.025,1);
    infec_Vero_pfu_95CI_max = quantile(outputvec_infec_Vero_pfu,0.975,1);
    infec_Vero_rna_95CI_min = quantile(outputvec_infec_Vero_rna,0.025,1);
    infec_Vero_rna_95CI_max = quantile(outputvec_infec_Vero_rna,0.975,1);


    CI95_min = {{decay_Vero_pfu_95CI_min;decay_Vero_rna_95CI_min};...
                    {infec_Vero_pfu_95CI_min;infec_Vero_rna_95CI_min}};
    CI95_max = {{decay_Vero_pfu_95CI_max;decay_Vero_rna_95CI_max};...
                    {infec_Vero_pfu_95CI_max;infec_Vero_rna_95CI_max}};

    %%% min and max %%%
    decay_pfu_min = min(outputvec_decay_Vero_pfu,[],1);
    decay_pfu_max = max(outputvec_decay_Vero_pfu,[],1);

    decay_rna_min = min(outputvec_decay_Vero_rna,[],1);
    decay_rna_max = max(outputvec_decay_Vero_rna,[],1);
    
    infec_pfu_min = min(outputvec_infec_Vero_pfu,[],1);
    infec_pfu_max = max(outputvec_infec_Vero_pfu,[],1);
    infec_rna_min = min(outputvec_infec_Vero_rna,[],1);
    infec_rna_max = max(outputvec_infec_Vero_rna,[],1);

    inf_min = {{decay_pfu_min;decay_rna_min};...
                {infec_pfu_min;infec_rna_min}};
    inf_max = {{decay_pfu_max;decay_rna_max};...
                {infec_pfu_max;infec_rna_max}};

    mlpars_Vero = readmatrix('maxLikValues.txt');
    ml_Vero_decay = decay_Vero(mlpars_Vero(2:end,2));
    ml_Vero_infec = model_Vero(mlpars_Vero(2:end,2));
    outputmaxlik = {ml_Vero_decay; ml_Vero_infec};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tiledplot = tiledlayout(1,2);
    set(gcf, 'Position',  [300, 100, 675, 305]);
    for aa = 1:2
        ax(aa) = nexttile(aa);
          set(ax(aa),...
          'box','on',...
          'XLim',XLims{aa},...
          'XTick',XTick{aa},...
          'YLim',YLims{aa},...   
          'YTick',YTicks{aa},...
          'FontSize',fontsize); 
            set(gca,'TickLength',[0.02, 0.01])
            xlabel('Time (hours post-infection)');
            ylabel('Log_{10} viral load per mL');
            hold on;
            for bb = 1:length(data{aa})
                patch(ax(aa),[t{aa},fliplr(t{aa})],[inf_min{aa}{bb},fliplr(inf_max{aa}{bb})],colors{aa,bb},'FaceAlpha',0.15,'EdgeColor',colors{aa,bb},'EdgeAlpha',0.15);
                patch(ax(aa),[t{aa},fliplr(t{aa})],[CI95_min{aa}{bb},fliplr(CI95_max{aa}{bb})],colors{aa,bb},'FaceAlpha',0.35,'EdgeColor',colors{aa,bb},'EdgeAlpha',0.35);
                hold on;
                plot(ax(aa),t{aa},outputmaxlik{aa}{bb},'Color',colors{aa,bb},'LineWidth',1,'LineStyle','-');
                hold on;
                data_var = {log10(data{aa}{bb}).Var2;
                            log10(data{aa}{bb}).Var3;
                            log10(data{aa}{bb}).Var4};
                data_std = stdev{aa}{bb}.std;
                mat_data_var = [transpose(data_var{1});transpose(data_var{2});transpose(data_var{3})];
                mean_mat_data_var = mean(mat_data_var,1);
                errorbar(time{aa},mean_mat_data_var,transpose(data_std),'o',...
                    'Color',colors{aa,bb},...
                    'MarkerSize',7,...
                    'MarkerEdgeColor',colors{aa,bb},...
                    'MarkerFaceColor','white',...
                    'CapSize',7,...
                    'LineWidth',1,...
                    'LineStyle','none');
            end
            h = [plot(NaN,NaN,'o','Color',colors{aa,1},'DisplayName','infectious (PFU)') 
                 plot(NaN,NaN,'o','Color',colors{aa,2},'DisplayName','total (RNA)')];
            set(h(1),'LineWidth',1);
            set(h(2),'LineWidth',1);
    
            legend(h,'Location','southwest');
    
            text(0.04,0.9,figlabs{aa},...
            'Units','Normalized',...
            'HorizontalAlignment','left',...
            'FontSize',fontsize,...
            'FontWeight','Normal');   
    end

    tiledplot.TileSpacing = 'compact';
    tiledplot.Padding = 'compact'; 

    exportgraphics(gcf,'./figures/kinetics_iol.png','Resolution',600);
    % exportgraphics(gcf,'../../LaTeX/figures/kinetics_iol.eps','ContentType','vector');


    %%% Helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%% Washing %%%
    function w = washing(t,w0)

        % w0 ... strength of washing
        % wd ... standard deviation of the length of washing (hours)
        % wt ... time of washing implementation (hours)
        wd = 0.05;
        wt = 1; 
        w = w0*1/sqrt(2*pi*wd^2)*exp(-(t-wt)^2/(2*wd^2));

    end

    %%% Decay %%%
    function [logV] = decay_Vero(par) 

        V_pfu = [];
        for ii = 1:par(13)
            V_pfu_ii = par(9)*(par(13)/par(6)*tt_decay).^(ii-1)/factorial(ii-1).*exp(-par(13)/par(6)*tt_decay);
            V_pfu = [V_pfu; V_pfu_ii];
        end
        logV_pfu = log10(sum(V_pfu,1));

        logV_rna = log10(par(10)*exp(-par(7)*tt_decay));

        logV = {logV_pfu;logV_rna};

    end
    
    %%% ODE system %%%
    function xdot = ode_Vero(t,x,par)

        %%% washing %%%
        omega = washing(t,par(8));

        %%% model %%%
        xdot = zeros(3+par(12)+par(13),1);

        xdot(1) = -par(1)*(sum(x(3+par(12):2+par(12)+par(13))))*x(1);    
        xdot(2) = par(1)*(sum(x(3+par(12):2+par(12)+par(13))))*x(1) - par(12)/par(2)*x(2);
        for ii = 3:1+par(12)
            xdot(ii) = par(12)/par(2)*(x(ii-1) - x(ii));
        end
        xdot(2+par(12)) = par(12)/par(2)*x(1+par(12)) - 1/par(3)*x(2+par(12));
        xdot(3+par(12)) = par(4)*x(2+par(12)) - omega*x(3+par(12)) - par(12)/par(6)*x(3+par(12));
        for kk = 4+par(12):2+par(12)+par(13)
           xdot(kk) = par(12)/par(6)*(x(kk-1) - x(kk)) - omega*x(kk);
        end      
        xdot(3+par(12)+par(13)) = par(5)*x(2+par(12)) - omega*x(3+par(12)+par(13)) - par(7)*x(3+par(12)+par(13));

    end

    %%% ode model solution %%%
    function [logV] = model_Vero(par) 

        par(12) = round(par(12));
        par(13) = round(par(13));

        sol = ode23s(@ode_Vero,[0 tmax],[T0 zeros(1,par(12)) 0 V0_pfu_Vero zeros(1,par(13)-1) par(11)],[],par);

        PFU = [];        
        for kk = 3+par(12):2+par(12)+par(13)
            V_kk = deval(sol,tt_infec,kk);
            PFU = [PFU; V_kk];
        end
        logPFU = log10(sum(PFU,1));
        logRNA = log10(deval(sol,tt_infec,3+par(12)+par(13)));

        logV={logPFU;logRNA};

    end
         

end