function [] = plot_data()

    close all
    clear all
    clc
    set(0,'DefaultFigureVisible','on');

    %%% generate an output directory
    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end

    %%% Vero cells CARRIB %%%
    data_CARRIB_PFU_Vero_low    = readtable('../../Supplementary_Data.xlsx','Sheet','CARRIB PFU Vero MOI 0.01');
    data_CARRIB_PFU_Vero_high   = readtable('../../Supplementary_Data.xlsx','Sheet','CARRIB PFU Vero MOI 1');
    data_CARRIB_RNA_Vero_low    = readtable('../../Supplementary_Data.xlsx','Sheet','CARRIB RNA Vero MOI 0.01');
    data_CARRIB_RNA_Vero_high   = readtable('../../Supplementary_Data.xlsx','Sheet','CARRIB RNA Vero MOI 1');
    data_CARRIB_PFU_Decay       = readtable('../../Supplementary_Data.xlsx','Sheet','CARRIB PFU Vero Decay');
    data_CARRIB_RNA_Decay       = readtable('../../Supplementary_Data.xlsx','Sheet','CARRIB RNA Vero Decay');

    data = {{data_CARRIB_PFU_Decay;
            data_CARRIB_RNA_Decay;};
            {data_CARRIB_PFU_Vero_low;
            data_CARRIB_RNA_Vero_low;};
            {data_CARRIB_PFU_Vero_high;
            data_CARRIB_RNA_Vero_high}};

    colors = {[0 0.5 0],[0 0 0];...
              [0 0 1],[0 0 0];...
              [1 0 0],[0 0 0]};
    YLims = {[-0.4 8.4];...
             [-0.5 10.5];...
             [-0.5 10.5]
             };   
    figlabs     = {'(A) Viral decay','(B) Infection 0.01 PFU/cell','(C) Infection 1 PFU/cell'};

    fontsize = 12;

    time_decay  = [0,8,24,48,72];
    time_infec  = [1,9,25,49,73];
    timeTick = {time_decay;...
                time_infec;...
                time_infec};
    xLims = {[-3 75];...
            [-2 76];...
            [-2 76]};


    tiledplot = tiledlayout(1,3,'TileSpacing','Compact');
    set(gcf, 'Position',  [300, 100, 1000, 300]);
    for aa = 1:length(data)
        ax(aa) = nexttile(aa);
        set(ax(aa),...
              'box','on',...
              'XLim',xLims{aa},...
              'XTick',timeTick{aa},...
              'YLim',YLims{aa},...
              'FontSize',fontsize); 
        set(gca,'TickLength',[0.02, 0.01])
        xlabel('Time (hours post-infection)');
        ylabel('Log_{10} viral load per mL');
        hold on;
        for bb = 1:length(data{aa})
            data_var = {log10(data{aa}{bb}.Var2);
                        log10(data{aa}{bb}.Var3);
                        log10(data{aa}{bb}.Var4)};
            for cc = 1:length(data_var)
                p=scatter(data{aa}{bb}.Var1,data_var{cc},50,'o',...
                    'MarkerEdgeColor',colors{aa,bb},...
                    'MarkerFaceColor','white',...
                    'LineWidth',1);
            end
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

    exportgraphics(gcf,'figures/data_carrib.png','Resolution',600);
    % exportgraphics(gcf,'../../LaTeX/figures/data_carrib.eps','ContentType','vector');

end