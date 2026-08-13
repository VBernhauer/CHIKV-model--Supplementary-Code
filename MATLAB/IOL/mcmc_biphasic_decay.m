function [] = mcmc_biphasic_decay()

    close all;
    clear all;
    clc;

    %%% 28 Degrees IOL Decay %%%
    data_IOL_PFU_Decay = readtable('../../Supplementary_Data.xlsx','Sheet','IOL PFU 28 Degrees Decay');
    data_IOL_RNA_Decay = readtable('../../Supplementary_Data.xlsx','Sheet','IOL RNA 28 Degrees Decay');

    stdev_IOL_PFU_Decay = mean(std(log10(data_IOL_PFU_Decay(:,2:end)),0,2,'omitnan'),'omitnan');
    stdev_IOL_RNA_Decay = mean(std(log10(data_IOL_RNA_Decay(:,2:end)),0,2,'omitnan'),'omitnan');

    data = {log10(data_IOL_PFU_Decay);
            log10(data_IOL_RNA_Decay)};
    stdev = {stdev_IOL_PFU_Decay;
             stdev_IOL_RNA_Decay};

    time_decay_plot  = [0,8,16,24,48,72,96,120];
    time_decay       = [0,2,4,8,16,24,48,72,96,120];
    tt_decay         = 0:0.1:120;

    c_pfu_high      = 2e-1;
    c_pfu           = 1e-2;
    c_pfu_low       = 1e-2;
    c_rna           = 2e-3;
    V0_pfu_stock    = 5e+5;
    V0_rna_stock    = 5e+7;
    par = log10([c_pfu_high,c_pfu,c_pfu_low,c_rna,V0_pfu_stock,V0_rna_stock]);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SLICE SAMPLER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nSamples    = 20000;
    burnIn      = 10000;
    thin        = 2;
    chains      = slicesample(par, nSamples, 'logpdf', @logprob);
    chains      = chains(burnIn+1:thin:end,:);

    %%% Bands %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    outputvec_decay_pfu = [];
    outputvec_decay_rna = [];
    for nn = 1:size(chains,1)
        
        if mod(nn,100) == 0
            disp(nn);
        end

        [logVdecay_out] = decay_tt(10.^chains(nn,:));

        outputvec_decay_pfu = [outputvec_decay_pfu; logVdecay_out{1}];
        outputvec_decay_rna = [outputvec_decay_rna; logVdecay_out{2}];  
        
    end

    %%% 95% CI %%%
    decay_Vero_pfu_95CI_min = quantile(outputvec_decay_pfu,0.025,1);
    decay_Vero_pfu_95CI_max = quantile(outputvec_decay_pfu,0.975,1);
    decay_Vero_rna_95CI_min = quantile(outputvec_decay_rna,0.025,1);
    decay_Vero_rna_95CI_max = quantile(outputvec_decay_rna,0.975,1);

    CI95_min = {decay_Vero_pfu_95CI_min;decay_Vero_rna_95CI_min};
    CI95_max = {decay_Vero_pfu_95CI_max;decay_Vero_rna_95CI_max};

    %%% min and max %%%
    decay_pfu_min = min(outputvec_decay_pfu,[],1);
    decay_pfu_max = max(outputvec_decay_pfu,[],1);

    decay_rna_min = min(outputvec_decay_rna,[],1);
    decay_rna_max = max(outputvec_decay_rna,[],1);   

    inf_min = {decay_pfu_min;decay_rna_min};
    inf_max = {decay_pfu_max;decay_rna_max};

    %%% ML value
    mlvalue = -99999.99;
    mlpars = zeros(1,size(chains,2));
    mlvalue_array = [];
    for ll = 1:size(chains,1)
        llvalue = loglike(10.^chains(ll,:));
        mlvalue_array = [mlvalue_array; llvalue];
        if llvalue > mlvalue
            mlvalue = llvalue;
            mlpars = chains(ll,:);
        end
    end
    outputmaxlik = decay_tt(10.^mlpars);

    %%% statistics %%%   
    a{1} = 'c_{pfu,high}';
    a{2} = 'c_{pfu}';
    a{3} = 'c_{pfu,low}';
    a{4} = 'c_{rna}';
    a{5} = 'V_{pfu,stock}';
    a{6} = 'V_{rna,stock}';
    labs = {a{1},...
            a{2},...
            a{3},...
            a{4},...
            a{5},...
            a{6}};

    names = {'parameter','mean','median','95CI lower', '95CI upper'};

    parsmean = 10.^mean(chains);
    parsmedian = 10.^median(chains);
    parsquantile = [];
    for nn = 1:size(chains,2)
        parsquantile = [parsquantile; quantile(10.^chains(:,nn),[0.025 0.975])];
    end  

    fid = fopen(strcat('posteriorValues_iol_28degrees.txt'),'w');
    fprintf(fid, '%2s %2s %2s %2s %2s\n', names{:});
    for nn = 1:length(labs)
        fprintf(fid,'%0s %.8f %.8f %.8f %.8f\n',labs{nn},parsmean(nn),parsmedian(nn),parsquantile(nn,1),parsquantile(nn,2));
    end
    fclose(fid);  

    %%% maximum likelihood values %%%
    fidml = fopen(strcat('maxLikValues_iol_28degrees.txt'),'w');
    fprintf(fidml, '%2s %.8f \n', 'MLvalue', mlvalue);
    for nn = 1:length(labs)
        fprintf(fidml,'%0s %.8f\n',labs{nn},mlpars(nn));
    end
    fclose(fidml);
    
    %%% Dynamics

    colors = {[0 0.5 0],[0 0 0]};
    YLims  = [2.75 8.25];   
    fontsize = 12;
    XLims = [-4 124];

    figure(1);
    set(gca,'box','on',...
          'XLim',XLims,...
          'XTick',time_decay_plot,...
          'YLim',YLims,...              
          'FontSize',fontsize); 
        set(gca,'TickLength',[0.02, 0.01])
        xlabel('Time (hours post-infection)');
        ylabel('Log_{10} viral load');
    hold on;
    for aa = 1:length(data)
        patch([tt_decay,fliplr(tt_decay)],[inf_min{aa},fliplr(inf_max{aa})],colors{aa},'FaceAlpha',0.15,'EdgeColor',colors{aa},'EdgeAlpha',0.15);
        hold on;
        patch([tt_decay,fliplr(tt_decay)],[CI95_min{aa},fliplr(CI95_max{aa})],colors{aa},'FaceAlpha',0.35,'EdgeColor',colors{aa},'EdgeAlpha',0.35);
        hold on;
        plot(tt_decay,outputmaxlik{aa},'Color',colors{aa},'LineWidth',1,'LineStyle','-');
        hold on;
    end

    for aa = 1:length(data)
        data_var = {data{aa}.Var2;
                    data{aa}.Var3;
                    data{aa}.Var4};
        for bb = 1:length(data_var)
            p=scatter(time_decay,data_var{bb},50,'o',...
                'MarkerEdgeColor',colors{aa},...
                'MarkerFaceColor','white');
        end
        hold on;
    end

    h = [plot(NaN,NaN,'o','Color',colors{1},'DisplayName','infectious (PFU/mL)') 
         plot(NaN,NaN,'o','Color',colors{2},'DisplayName','total (RNA/mL)')];
    legend(h,'Location','southwest');

    %%% save figure %%%
    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end

    exportgraphics(gcf,'./figures/kinetics_iol_28degrees.png','Resolution',600);
    % exportgraphics(gcf,'../../LaTeX/figures/kinetics_iol_28degrees.eps','ContentType','vector');


    %%% Posterior histograms

    a{1} = 'Log_{10} c_{pfu}^{high}';
    a{2} = 'Log_{10} c_{pfu}';
    a{3} = 'Log_{10} c_{pfu}^{low}';
    a{4} = 'Log_{10} c_{rna}';
    a{5} = 'Log_{10} V_{pfu}^{stock}';
    a{6} = 'Log_{10} V_{rna}^{stock}';
    figlabs = {a{1},...
            a{2},...
            a{3},...
            a{4},...
            a{5},...
            a{6}};

    fontsize = 12;

    XLim = {[-2 -0.5],...
            [-1.5 0],...
            [-2.4 -1.4],...
            [-3.0 -2.2],...
            [5.4 6.0],...
            [7.60 7.75]};
    XTicks = {[-2,-1.5,-1,-0.5],...
              [-1.5,-1,-0.5,0],...
              [-2.4,-1.9,-1.4],...
              [-3.0,-2.6,-2.2],...
              [5.4,5.6,5.8,6.0],...
              [7.6,7.65,7.7,7.75]};
              
    figure(2);
    tiledplot = tiledlayout(2,3,'TileSpacing','Compact');
    set(gcf, 'Position',  [300, 100, 800, 400]);
    for aa = 1:length(figlabs)
        ax(aa) = nexttile(aa);
        set(ax(aa),...
              'box','on',...
              'XLim', XLim{aa},...
              'XTick', XTicks{aa},...
              'YLim', [0 0.2],...
              'FontSize',fontsize); 
        set(gca,'TickLength',[0.02, 0.01])
        hold on;
        histogram(chains(:,aa), Normalization="probability",...
                                DisplayStyle="bar",...
                                FaceColor=[0 0 0],...
                                FaceAlpha=0.35,...
                                EdgeColor=[0 0 0],...
                                EdgeAlpha=1);   
        text(0.075,0.85,figlabs{aa},...
            'Units','Normalized',...
            'HorizontalAlignment','left',...
            'FontSize',fontsize,...
            'FontWeight','Normal');
    end
    ylabel(tiledplot,'Normalized frequency', fontsize=fontsize);
    tiledplot.TileSpacing = 'compact';
    tiledplot.Padding = 'compact'; 

    exportgraphics(gcf,'./figures/posteriors_iol_28degrees.png','Resolution',600);
    % exportgraphics(gcf,'../LaTeX/figures/posteriors_iol_28degrees.eps','ContentType','vector');
   

    %%% Helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%% Decay %%%
    function [logV] = decay_tt(par) 

        logV_pfu = log10(par(5)*exp(-(par(1)+par(2))*tt_decay).*(1-par(1)/(par(1)+par(2)-par(3))*(1-exp((par(1)+par(2)-par(3))*tt_decay))));
        % logV_pfu = log10(par(5)*exp(-par(1)*tt_decay));
        logV_rna = log10(par(6)*exp(-par(4)*tt_decay));

        logV = {logV_pfu;logV_rna};

        end

    function [logV] = decay(par) 

        logV_pfu = log10(par(5)*exp(-(par(1)+par(2))*time_decay).*(1-par(1)/(par(1)+par(2)-par(3))*(1-exp((par(1)+par(2)-par(3))*time_decay))));
        % logV_pfu = log10(par(5)*exp(-par(1)*time_decay));
        logV_rna = log10(par(6)*exp(-par(4)*time_decay));

        logV = {logV_pfu;logV_rna};

    end
         
    %%% log-likelihood %%%
    function [value] = loglike(par)  

        logV  = decay(par);
        
        value = 0;
        for ii = 1:length(logV)      
            logdata  = [data{ii}.Var2,...
                        data{ii}.Var3,...
                        data{ii}.Var4];
            logstdev = stdev{ii}.std;
            logmodel = logV{ii}.';
            for kk = 1:3
                idx=find(~isnan(logdata(:,kk)));
                value = value+sum(lognormpdf(logdata(idx,kk),logmodel(idx,1),logstdev));
            end
        end

    end

    %%% log-likelihood formula %%%
    function [value] = lognormpdf(logdata,logmodel,logstdev)
        value = -0.5*((logdata-logmodel)./logstdev).^2  - log(sqrt(2*pi).*logstdev); 
    end

    %%% priors %%%
    function [flag] = logprior(par)
        if all(par>0) && par(1)>par(3)
            flag = 0;
        else
            flag = -inf;
        end       
    end

    %%% log-probability %%%
    function [value] = logprob(par)
        par = 10.^par;
        lp = logprior(par);
        if ~isfinite(lp)
            value = -inf;
        else
            value = lp+loglike(par);
        end       
    end

   
end