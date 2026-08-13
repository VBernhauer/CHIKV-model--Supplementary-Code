function [] = plot_practical_identifiability()

    clc;
    close all;
    MLvalues        = load(strcat('./MLvalues.mat'));

    %%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    a{1} = 'Log_{10} \beta';
    a{2} = 'Log_{10} c_{rna}';
    parlabs = {a{1},...
               a{2}};

    figlabs = {'(B)','(C)'};

    parameters = {'beta',...
                  'c_rna'};
    idx = [1,7];

    XLabs = {[-5,-4,-3,-2,-1,0];...
             [-3,-2.5,-2]};
    YLims = {[-10 40];...
             [-10 15]};

    fontsize = 12;

    tiledplot = tiledlayout(2,1,'TileSpacing','Compact');
    set(gcf, 'Position',  [300, 100, 250, 500]);
    for aa = 1:length(parlabs)
        fval = load(strcat('./PIA/',parameters{aa},'.mat'));
        xXaxis = log10(fval.funval(:,1));
        yYaxis = fval.funval(:,2);
        ax(aa) = nexttile(aa);  
        set(ax(aa),...
          'box','on',...
          'XLim',[min(xXaxis), max(xXaxis)],...
          'XTick',XLabs{aa},...
          'YLim',YLims{aa},...
          'FontSize',fontsize); 
        set(gca,'TickLength',[0.02, 0.01])
        hold on;
        plot(ax(aa),xXaxis,yYaxis,...
                                'LineWidth',1,...
                                'Color',[0 0 0],...
                                'LineStyle','-')
        hold on
        plot(ax(aa),xXaxis,(-MLvalues.MLvalues(1,1)+0.5*chi2inv(0.95,1))*ones(1,length(xXaxis)),'Color',[0 0 0 1],'LineWidth',1,'LineStyle','--')
        hold on
        plot(ax(aa),log10(MLvalues.MLvalues(1,idx(aa)+1)),-MLvalues.MLvalues(1,1),'Marker','o',...
                    'MarkerSize',8,...
                    'MarkerFaceColor',[1 0 0],...
                    'MarkerEdgeColor',[1 0 0],...
                    'Color',[1 0 0],...
                    'LineStyle','none')
        xlabel(ax(aa),parlabs{aa},'FontSize',fontsize);
        text(0.05,0.9,figlabs{aa},...
            'Units','Normalized',...
            'HorizontalAlignment','left',...
            'FontSize',fontsize,...
            'FontWeight','Normal');
    hold on;
    end    
    ylabel(tiledplot,'Negative profile likelihood','FontSize',fontsize);
        
    tiledplot.TileSpacing = 'compact';
    tiledplot.Padding = 'compact';

    exportgraphics(gcf,strcat('../figures/practical_identifiability_beta_crna.png'));
    exportgraphics(gcf,strcat('../../../LaTeX/figures/practical_identifiability_beta_crna.eps'),'ContentType','vector');


end