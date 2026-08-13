function [] = mcmc_ml_95ci()

    %%% command:
    %%% mcmc_ml_95ci(5)
   
    %%% Vero cells IOL %%%
    data_IOL_PFU_Vero       = readtable('../../Supplementary_Data.xlsx','Sheet','IOL PFU Vero');
    data_IOL_RNA_Vero       = readtable('../../Supplementary_Data.xlsx','Sheet','IOL RNA Vero');
    data_IOL_PFU_Vero_Decay = readtable('../../Supplementary_Data.xlsx','Sheet','IOL PFU Vero Decay');
    data_IOL_RNA_Vero_Decay = readtable('../../Supplementary_Data.xlsx','Sheet','IOL RNA Vero Decay');
    
    stdev_IOL_PFU_Vero       = mean(std(log10(data_IOL_PFU_Vero(:,2:end)),0,2,'omitnan'),'omitnan');
    stdev_IOL_RNA_Vero       = mean(std(log10(data_IOL_RNA_Vero(:,2:end)),0,2,'omitnan'),'omitnan');
    stdev_IOL_PFU_Vero_Decay = mean(std(log10(data_IOL_PFU_Vero_Decay(:,2:end)),0,2,'omitnan'),'omitnan');
    stdev_IOL_RNA_Vero_Decay = mean(std(log10(data_IOL_RNA_Vero_Decay(:,2:end)),0,2,'omitnan'),'omitnan');

    tmax        = 121;
    time_decay  = [0,2,4,8,16,24,48,72,96,120];
    time        = [1.5,3,5,9,17,25,49,73,97,121];

    %%% global fixed parameters
    T0      = 1;    

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

    chains = [];
    for nn = 1:5
        kkchains = load(strcat('./chains/chains_',num2str(nn),'.mat'));
        chains = [chains; kkchains.chains(:,:);];
    end

    mlvalue = -99999.99;
    mlpars = zeros(1,size(chains,2));
    mlvalue_array = [];
    for ll = 1:size(chains,1)
        if mod(ll,100) == 0
            disp(ll);
        end
        llvalue = loglike(10.^chains(ll,:));
        mlvalue_array = [mlvalue_array; llvalue];
        if llvalue > mlvalue
            mlvalue = llvalue;
            mlpars = 10.^chains(ll,:);
        end
    end

    a{1} = 'beta';
    a{2} = 'tau_L';
    a{3} = 'tau_I';
    a{4} = 'p_pfu';
    a{5} = 'p_rna';
    a{6} = 'tau_V';
    a{7} = 'c_rna';
    a{8} = 'omega_0';
    a{9} = 'V0_pfu_stock';
    a{10} = 'V0_rna_stock';
    a{11} = 'V0_rna';
    a{12} = 'n_L';
    a{13} = 'n_V';
    labs = {a{1},...
            a{2},...
            a{3},...
            a{4},...
            a{5},...
            a{6},...
            a{7},...
            a{8},...
            a{9},...
            a{10},...
            a{11},...
            a{12},...
            a{13}};        

    %%% save chains %%%
    save(strcat('./mlvalues.mat'),'mlvalue_array');
            
    %%% statistics %%%    
    names = {'parameter','mean','median','95CI lower', '95CI upper'};

    parsmean = 10.^mean(chains);
    parsmedian = 10.^median(chains);
    parsquantile = [];
    for nn = 1:size(chains,2)
        parsquantile = [parsquantile; quantile(10.^chains(:,nn),[0.025 0.975])];
    end  

    fid = fopen(strcat('posteriorValues.txt'),'w');
    fprintf(fid, '%2s %2s %2s %2s %2s\n', names{:});
    for nn = 1:length(labs)
        fprintf(fid,'%0s %.8f %.8f %.8f %.8f\n',labs{nn},parsmean(nn),parsmedian(nn),parsquantile(nn,1),parsquantile(nn,2));
    end
    fclose(fid);  

    %%% maximum likelihood values %%%
    fidml = fopen(strcat('maxLikValues.txt'),'w');
    fprintf(fidml, '%2s %.8f \n', 'MLvalue', mlvalue);
    mlpars(end-1) = round(mlpars(end-1));
    mlpars(end)   = round(mlpars(end));
    for nn = 1:length(labs)
        fprintf(fidml,'%0s %.8f\n',labs{nn},mlpars(nn));
    end
    fclose(fidml);
    
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
            V_pfu_ii = par(9)*(par(13)/par(6)*time_decay).^(ii-1)/factorial(ii-1).*exp(-par(13)/par(6)*time_decay);
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
            V_kk = deval(sol,time,kk);
            PFU = [PFU; V_kk];
        end
        logPFU = log10(sum(PFU,1));
        logRNA = log10(deval(sol,time,3+par(12)+par(13)));

        logV={logPFU;logRNA};

    end
         
    %%% Log-likelihood %%%
    function [value] = loglike(par)  

        logV_decay  = decay_Vero(par);
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

    %%% Log-likelihood formula %%%
    function [value] = lognormpdf(logdata,logmodel,logstdev)
        value = -0.5*((logdata-logmodel)./logstdev).^2  - log(sqrt(2*pi).*logstdev); 
    end

end