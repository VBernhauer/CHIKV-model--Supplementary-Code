function [] = plot_mcmc()

    close all
    clear all
    clc
    set(0,'DefaultFigureVisible','on');

    %%% generate an output directory
    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end

    %%% Vero cells CARRIB data %%%
    data_CARRIB_PFU_Decay       = readtable('../../Supplementary_Data.xlsx','Sheet','CARRIB PFU Vero Decay');
    data_CARRIB_RNA_Decay       = readtable('../../Supplementary_Data.xlsx','Sheet','CARRIB RNA Vero Decay');
    data_CARRIB_PFU_Vero_low    = readtable('../../Supplementary_Data.xlsx','Sheet','CARRIB PFU Vero MOI 0.01');
    data_CARRIB_PFU_Vero_high   = readtable('../../Supplementary_Data.xlsx','Sheet','CARRIB PFU Vero MOI 1');
    data_CARRIB_RNA_Vero_low    = readtable('../../Supplementary_Data.xlsx','Sheet','CARRIB RNA Vero MOI 0.01');
    data_CARRIB_RNA_Vero_high   = readtable('../../Supplementary_Data.xlsx','Sheet','CARRIB RNA Vero MOI 1');

    data = {{log10(data_CARRIB_PFU_Decay);
            log10(data_CARRIB_RNA_Decay)};
            {log10(data_CARRIB_PFU_Vero_low);
            log10(data_CARRIB_RNA_Vero_low)};
            {log10(data_CARRIB_PFU_Vero_high);
            log10(data_CARRIB_RNA_Vero_high)}};

    stdev_CARRIB_PFU_Decay      = std(log10(data_CARRIB_PFU_Decay(:,2:end)),0,2,'omitnan');
    stdev_CARRIB_RNA_Decay      = std(log10(data_CARRIB_RNA_Decay(:,2:end)),0,2,'omitnan');
    stdev_CARRIB_PFU_Vero_low   = std(log10(data_CARRIB_PFU_Vero_low(:,2:end)),0,2,'omitnan');
    stdev_CARRIB_PFU_Vero_high  = std(log10(data_CARRIB_PFU_Vero_high(:,2:end)),0,2,'omitnan');
    stdev_CARRIB_RNA_Vero_low   = std(log10(data_CARRIB_RNA_Vero_low(:,2:end)),0,2,'omitnan');
    stdev_CARRIB_RNA_Vero_high  = std(log10(data_CARRIB_RNA_Vero_high(:,2:end)),0,2,'omitnan');

    stdev = {{stdev_CARRIB_PFU_Decay;
             stdev_CARRIB_RNA_Decay};
             {stdev_CARRIB_PFU_Vero_low;
             stdev_CARRIB_PFU_Vero_high};
             {stdev_CARRIB_RNA_Vero_low;
             stdev_CARRIB_RNA_Vero_high}};

    colors = {[0 0.5 0],[0 0 0];...
              [0 0 1],[0 0 0];...
              [1 0 0],[0 0 0]};
    YLims = {[-0.4 8.4];...
             [-0.5 10.5];...
             [-0.5 10.5]
             }; 
    XTick = {[0,8,24,48,72];...
             [1,9,25,49,73];...
             [1,9,25,49,73]};
    XLims = {[-3 75];...
             [-2 76];...
             [-2 76]};
    figlabs     = {'(A) Viral decay','(B) Infection 0.01 PFU/cell','(C) Infection 1 PFU/cell'};
    fontsize = 12;

    tmax    = 73;
    T0      = 1;

    time_decay  = [0,4,6,8,24,48,72];
    time_infec  = [1.5,5,7,9,25,49,73];
    time = {time_decay;...
            time_infec;...
            time_infec};

    tt_decay = 0:0.1:72;
    tt_infec = 0:0.1:73;
    t = {tt_decay;...
         tt_infec;...
         tt_infec};
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    outputvec_decay_pfu         = readmatrix('output_mcmc_CARRIB_decay_pfu.txt');
    outputvec_decay_rna         = readmatrix('output_mcmc_CARRIB_decay_rna.txt');
    outputvec_infec_pfu_low     = readmatrix('output_mcmc_CARRIB_pfu_001.txt');
    outputvec_infec_rna_low     = readmatrix('output_mcmc_CARRIB_rna_001.txt');
    outputvec_infec_pfu_high    = readmatrix('output_mcmc_CARRIB_pfu_1.txt');
    outputvec_infec_rna_high    = readmatrix('output_mcmc_CARRIB_rna_1.txt');

    decay_pfu_95CI_min = quantile(outputvec_decay_pfu,0.025,1);
    decay_pfu_95CI_max = quantile(outputvec_decay_pfu,0.975,1);

    decay_rna_95CI_min = quantile(outputvec_decay_rna,0.025,1);
    decay_rna_95CI_max = quantile(outputvec_decay_rna,0.975,1);

    infec_pfu_low_95CI_min = quantile(outputvec_infec_pfu_low,0.025,1);
    infec_pfu_low_95CI_max = quantile(outputvec_infec_pfu_low,0.975,1);
    infec_rna_low_95CI_min = quantile(outputvec_infec_rna_low,0.025,1);
    infec_rna_low_95CI_max = quantile(outputvec_infec_rna_low,0.975,1);

    infec_pfu_high_95CI_min = quantile(outputvec_infec_pfu_high,0.025,1);
    infec_pfu_high_95CI_max = quantile(outputvec_infec_pfu_high,0.975,1);
    infec_rna_high_95CI_min = quantile(outputvec_infec_rna_high,0.025,1);
    infec_rna_high_95CI_max = quantile(outputvec_infec_rna_high,0.975,1);

    CI95_min = {{decay_pfu_95CI_min;decay_rna_95CI_min};...
                {infec_pfu_low_95CI_min;infec_rna_low_95CI_min};...
                {infec_pfu_high_95CI_min;infec_rna_high_95CI_min}};
    CI95_max = {{decay_pfu_95CI_max;decay_rna_95CI_max};...
                {infec_pfu_low_95CI_max;infec_rna_low_95CI_max};...
                {infec_pfu_high_95CI_max;infec_rna_high_95CI_max}};

    decay_pfu_min = min(outputvec_decay_pfu,[],1);
    decay_pfu_max = max(outputvec_decay_pfu,[],1);

    decay_rna_min = min(outputvec_decay_rna,[],1);
    decay_rna_max = max(outputvec_decay_rna,[],1);

    infec_pfu_low_min = min(outputvec_infec_pfu_low,[],1);
    infec_pfu_low_max = max(outputvec_infec_pfu_low,[],1);
    infec_rna_low_min = min(outputvec_infec_rna_low,[],1);
    infec_rna_low_max = max(outputvec_infec_rna_low,[],1);

    infec_pfu_high_min = min(outputvec_infec_pfu_high,[],1);
    infec_pfu_high_max = max(outputvec_infec_pfu_high,[],1);
    infec_rna_high_min = min(outputvec_infec_rna_high,[],1);
    infec_rna_high_max = max(outputvec_infec_rna_high,[],1);

    inf_min = {{decay_pfu_min;decay_rna_min};...
                {infec_pfu_low_min;infec_rna_low_min};...
                {infec_pfu_high_min;infec_rna_high_min}};
    inf_max = {{decay_pfu_max;decay_rna_max};...
                {infec_pfu_low_max;infec_rna_low_max};...
                {infec_pfu_high_max;infec_rna_high_max}};

    mlpars = readmatrix('maxLikValues.txt');
    ml_decay = decay(10.^mlpars(2:end,2));
    ml_infec = model(10.^mlpars(2:end,2));
    outputmaxlik = {ml_decay; ml_infec{1}; ml_infec{2}};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tiledplot = tiledlayout(1,3,'TileSpacing','Compact');
    set(gcf, 'Position',  [300, 100, 1000, 300]);
    for aa = 1:length(data)
        ax(aa) = nexttile(aa);
        set(ax(aa),...
              'box','on',...
              'XLim',XLims{aa},...
              'XTick',XTick{aa},...
              'YLim',YLims{aa},...
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
            data_var = {data{aa}{bb}.Var2;
                        data{aa}{bb}.Var3;
                        data{aa}{bb}.Var4};
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

        legend(h,'Location','southeast');
        text(0.04,0.9,figlabs{aa},...
        'Units','Normalized',...
        'HorizontalAlignment','left',...
        'FontSize',fontsize,...
        'FontWeight','Normal');   
    end

    tiledplot.TileSpacing = 'compact';
    tiledplot.Padding = 'compact'; 

    exportgraphics(gcf,'figures/kinetics_carrib.png','Resolution',600);
    % exportgraphics(gcf,'../../LaTeX/figures/kinetics_carrib.eps','ContentType','vector');


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
    function [logV] = decay(par) 

        logV_pfu = log10(par(8)*exp(-par(6)*tt_decay));
        logV_rna = log10(par(9)*ones(1,length(tt_decay)));

        logV = {logV_pfu;logV_rna};

    end
    
    %%% ODE system %%%
    function xdot = ode(t,x,par)

        %%% washing %%%
        omega = washing(t,par(7));

        %%% model %%%
        xdot = zeros(4+par(10),1);

        xdot(1) = -par(1)*x(3+par(10))*x(1);    
        xdot(2) = par(1)*x(3+par(10))*x(1) - par(10)/par(2)*x(2);
        for ii = 3:1+par(10)
            xdot(ii) = par(10)/par(2)*(x(ii-1) - x(ii));
        end
        xdot(2+par(10)) = par(10)/par(2)*x(1+par(10)) - 1/par(3)*x(2+par(10));
        xdot(3+par(10)) = par(4)*x(2+par(10)) - omega*x(3+par(10)) - par(6)*x(3+par(10));
        xdot(4+par(10)) = par(5)*x(2+par(10)) - omega*x(4+par(10));

    end

    %%% ode model solution %%%
    function [logV] = model(par) 

        par(10) = round(par(10));
        % par(11) = round(par(11));

        %%% low MOI        
        sol_low = ode23s(@ode,[0 tmax],[T0 zeros(1,par(10)) 0 0.01*par(8) 0.01*par(9)],[],par);
        logPFU_low = log10(deval(sol_low,tt_infec,3+par(10)));
        logRNA_low = log10(deval(sol_low,tt_infec,4+par(10)));

        %%% high MOI
        sol_high = ode23s(@ode,[0 tmax],[T0 zeros(1,par(10)) 0 par(8) par(9)],[],par);
        logPFU_high = log10(deval(sol_high,tt_infec,3+par(10)));
        logRNA_high = log10(deval(sol_high,tt_infec,4+par(10)));

        logV={{logPFU_low;logRNA_low},{logPFU_high;logRNA_high}};

    end
         

end