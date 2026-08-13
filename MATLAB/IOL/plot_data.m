function [] = plot_data()

    close all
    clear all
    clc
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
    figlabs     = {'(A) Viral decay','(B) Infection 0.1 PFU/cell'};

    fontsize = 12;

    time_decay  = [0,8,24,48,72,96,120];
    time_infec  = [1,9,25,49,73,97,121];
    timeTick = {time_decay;...
                time_infec};
    XLims = {[-4 124];...
             [-3 125]};


    tiledplot = tiledlayout(1,2,'TileSpacing','Compact');
    set(gcf, 'Position',  [300, 100, 675, 305]);
    for aa = 1:length(figlabs)
        ax(aa) = nexttile(aa);
        set(ax(aa),...
              'box','on',...
              'XLim',XLims{aa},...
              'XTick',timeTick{aa},...
              'YLim',YLims{aa},... 
              'YTick',YTicks{aa},...
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
        legend(h,'Location','southwest');
        
        text(0.05,0.9,figlabs{aa},...
        'Units','Normalized',...
        'HorizontalAlignment','left',...
        'FontSize',fontsize,...
        'FontWeight','Normal');   
    end

    tiledplot.TileSpacing = 'compact';
    tiledplot.Padding = 'compact'; 

    exportgraphics(gcf,'./figures/data_iol.png','Resolution',600);
    % exportgraphics(gcf,'../../LaTeX/figures/data_iol.eps','ContentType','vector');

end