function [results] = mcmc(nstart,nend,par)

    %%% Vero cells CARRIB data %%%
    data_CARRIB_PFU_Decay       = readtable('../../Supplementary_Data.xlsx','Sheet','CARRIB PFU Vero Decay');
    data_CARRIB_RNA_Decay       = readtable('../../Supplementary_Data.xlsx','Sheet','CARRIB RNA Vero Decay');
    data_CARRIB_PFU_Vero_low    = readtable('../../Supplementary_Data.xlsx','Sheet','CARRIB PFU Vero MOI 0.01');
    data_CARRIB_PFU_Vero_high   = readtable('../../Supplementary_Data.xlsx','Sheet','CARRIB PFU Vero MOI 1');
    data_CARRIB_RNA_Vero_low    = readtable('../../Supplementary_Data.xlsx','Sheet','CARRIB RNA Vero MOI 0.01');
    data_CARRIB_RNA_Vero_high   = readtable('../../Supplementary_Data.xlsx','Sheet','CARRIB RNA Vero MOI 1');

    %%% Vero cells CARRIB standard deviation %%%
    stdev_CARRIB_PFU_Decay      = mean(std(log10(data_CARRIB_PFU_Decay(:,2:end)),0,2,'omitnan'),'omitnan');
    stdev_CARRIB_RNA_Decay      = mean(std(log10(data_CARRIB_RNA_Decay(:,2:end)),0,2,'omitnan'),'omitnan');
    stdev_CARRIB_PFU_Vero_low   = mean(std(log10(data_CARRIB_PFU_Vero_low(:,2:end)),0,2,'omitnan'),'omitnan');
    stdev_CARRIB_PFU_Vero_high  = mean(std(log10(data_CARRIB_PFU_Vero_high(:,2:end)),0,2,'omitnan'),'omitnan');
    stdev_CARRIB_RNA_Vero_low   = mean(std(log10(data_CARRIB_RNA_Vero_low(:,2:end)),0,2,'omitnan'),'omitnan');
    stdev_CARRIB_RNA_Vero_high  = mean(std(log10(data_CARRIB_RNA_Vero_high(:,2:end)),0,2,'omitnan'),'omitnan');

    data_decay = {log10(data_CARRIB_PFU_Decay);
                  log10(data_CARRIB_RNA_Decay)};

    stdev_decay = {stdev_CARRIB_PFU_Decay;
                   stdev_CARRIB_RNA_Decay};

    data = {{log10(data_CARRIB_PFU_Vero_low);
            log10(data_CARRIB_RNA_Vero_low)};
            {log10(data_CARRIB_PFU_Vero_high);
            log10(data_CARRIB_RNA_Vero_high)}};

    stdev = {{stdev_CARRIB_PFU_Vero_low;
            stdev_CARRIB_PFU_Vero_high};
            {stdev_CARRIB_RNA_Vero_low;
            stdev_CARRIB_RNA_Vero_high}};

    tmax        = 73;
    time_decay  = [0,4,6,8,24,48,72];
    time        = [1.5,5,7,9,25,49,73];


    %%% global fixed parameters
    T0      = 1;
    nEmin   = 1;
    nEmax   = 101;
    omegamax = 10;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SLICE SAMPLER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iround=nstart:nend

        nSamples    = 20000;
        burnIn      = 10000;
        thin        = 10;
        chains          = slicesample(par(iround,:), nSamples, 'logpdf', @logprob);
        chains          = chains(burnIn+1:thin:end,:);
        results{iround} = chains;

        %%% save chains %%%
        save(strcat('./chains/chains_',num2str(iround),'.mat'),'chains');

    end

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

        logV_pfu = log10(par(8)*exp(-par(6)*time_decay));
        logV_rna = log10(par(9)*ones(1,length(time_decay)));

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
        sol_low = ode23s(@ode,[0 tmax],[T0 zeros(1,par(10)) 0 0.01*par(8) 0.01*par(9)],[],par);
        try
            logPFU_low = log10(deval(sol_low,time,3+par(10)));
        catch
            logPFU_low = ones(1,length(time));
        end

        try
            logRNA_low = log10(deval(sol_low,time,4+par(10)));
        catch
            logRNA_low = ones(1,length(time));        
        end
        

        %%% high MOI
        sol_high = ode23s(@ode,[0 tmax],[T0 zeros(1,par(10)) 0 par(8) par(9)],[],par);
        try
            logPFU_high = log10(deval(sol_high,time,3+par(10)));
        catch
            logPFU_high = ones(1,length(time));
        end

        try
            logRNA_high = log10(deval(sol_high,time,4+par(10)));
        catch
            logRNA_high = ones(1,length(time));
        end

        logV={{logPFU_low;logRNA_low},{logPFU_high;logRNA_high}};

    end
         
    %%% log-likelihood %%%
    function [value] = loglike(par)  

        [logV] = model(par);
        value = 0;
        for ii = 1:length(logV)      
            for jj = 1:length(logV{ii})
                logdata  = [data{ii}{jj}.Var2,...
                            data{ii}{jj}.Var3,...
                            data{ii}{jj}.Var4];
                logstdev = stdev{ii}{jj}.std;
                logmodel = logV{ii}{jj}.';
                for kk = 1:3
                    idx=find(~isnan(logdata(:,kk)));
                    value = value+sum(lognormpdf(logdata(idx,kk),logmodel(idx,1),logstdev));
                end
            end
        end

        [logVdecay] = decay(par);
        for ii = 1:length(logVdecay) 
            logdata  = [data_decay{ii}.Var2,...
                        data_decay{ii}.Var3,...
                        data_decay{ii}.Var4];
            logstdev = stdev_decay{ii}.std;
            logmodel = logVdecay{ii}.';
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
        if all(par>0) && par(10)>=nEmin && par(10)<nEmax && par(7)<omegamax
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