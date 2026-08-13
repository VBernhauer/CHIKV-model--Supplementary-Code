function  [results] = mcmc(nstart,nend,par)

    %%% Vero cells IOL %%%
    data_IOL_PFU_Vero       = readtable('../../Supplementary_Data.xlsx','Sheet','IOL PFU Vero');
    data_IOL_RNA_Vero       = readtable('../../Supplementary_Data.xlsx','Sheet','IOL RNA Vero');
    data_IOL_PFU_Vero_Decay = readtable('../../Supplementary_Data.xlsx','Sheet','IOL PFU Vero Decay');
    data_IOL_RNA_Vero_Decay = readtable('../../Supplementary_Data.xlsx','Sheet','IOL RNA Vero Decay');
    
    stdev_IOL_PFU_Vero       = mean(std(log10(data_IOL_PFU_Vero(:,2:end)),0,2,'omitnan'),'omitnan');
    stdev_IOL_RNA_Vero       = mean(std(log10(data_IOL_RNA_Vero(:,2:end)),0,2,'omitnan'),'omitnan');  
    stdev_IOL_PFU_Vero_Decay = mean(std(log10(data_IOL_PFU_Vero_Decay(:,2:end)),0,2,'omitnan'),'omitnan');  
    stdev_IOL_RNA_Vero_Decay = mean(std(log10(data_IOL_RNA_Vero_Decay(:,2:end)),0,2,'omitnan'),'omitnan');

    data = {{log10(data_IOL_PFU_Vero_Decay);
            log10(data_IOL_RNA_Vero_Decay)};
            {log10(data_IOL_PFU_Vero);
            log10(data_IOL_RNA_Vero)}};

    stdev = {{stdev_IOL_PFU_Vero_Decay;
            stdev_IOL_RNA_Vero_Decay};
            {stdev_IOL_PFU_Vero;
            stdev_IOL_RNA_Vero}};

    T0_pfu  = 3e+5;
    MOI     = 0.1;
    V0_pfu  = MOI*T0_pfu;

    tmax        = 121;
    time_decay  = [0,2,4,8,16,24,48,72,96,120];
    time        = [1.5,3,5,9,17,25,49,73,97,121];


    %%% global fixed parameters
    T0          = 1;
    nEmin       = 1;
    nEmax       = 51;
    nVmin       = 1;
    nVmax       = 51;
    crnamax     = 5e-4;
    omegamax    = 10;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SLICE SAMPLER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nSamples    = 20000;
    burnIn      = 10000;
    thin        = 10;
    for iround=nstart:nend
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

        V_pfu = [];
        for ii = 1:par(13)
            V_pfu_ii = par(9)*exp(-par(13)/par(6)*time_decay).*((par(13)/par(6)*time_decay).^(ii-1))/factorial(ii-1);
            V_pfu = [V_pfu; V_pfu_ii];
        end
        logV_pfu = log10(sum(V_pfu,1));

        logV_rna = log10(par(10)*exp(-par(7)*time_decay));

        logV = {logV_pfu;logV_rna};

    end
    
    %%% ODE system %%%
    function xdot = ode(t,x,par)

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

    %%% ODE model solution %%%
    function [logV] = model(par) 

        par(12) = round(par(12));
        par(13) = round(par(13));

        sol = ode23s(@ode,[0 tmax],[T0 zeros(1,par(12)) 0 V0_pfu zeros(1,par(13)-1) par(11)],[],par);
        
        PFU = [];        
        for kk = 3+par(12):2+par(12)+par(13)
            try
                V_kk = deval(sol,time,kk);
            catch
                V_kk = 10*ones(1,length(time));
            end
            PFU = [PFU; V_kk];
        end
        logPFU = log10(sum(PFU,1));

        try
            logRNA = log10(deval(sol,time,3+par(12)+par(13)));
        catch
            logRNA = ones(1,length(time));        
        end

        logV={logPFU;logRNA};

    end
         
    %%% Log-likelihood %%%
    function [value] = loglike(par)  

        logV_decay  = decay(par);
        logV_inf    = model(par);
        logV = {logV_decay;logV_inf};
        
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

    end

    %%% Priors %%%
    function [flag] = logprior(par)
        if all(par>0) && par(7)>=crnamax ...
                      && par(8)<omegamax ...
                      && par(12)>=nEmin && par(12)<nEmax ...
                      && par(13)>=nVmin && par(13)<nVmax 
            flag = 0;
        else
            flag = -inf;
        end       
    end

    %%% Log-probability %%%
    function [value] = logprob(par)
        par = 10.^par;
        lp = logprior(par);
        if ~isfinite(lp)
            value = -inf;
        else
            value = lp+loglike(par);
        end       
    end

    %%% Log-likelihood formula %%%
    function [value] = lognormpdf(logdata,logmodel,logstdev)
        value = -0.5*((logdata-logmodel)./logstdev).^2  - log(sqrt(2*pi).*logstdev); 
    end

   
end