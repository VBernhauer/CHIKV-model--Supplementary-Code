function []=mcmc_params(n)
    %%% command:
    %%% mcmc_params(n)
    
    beta            = [1e-6 1e-5];%1
    tau_L           = [1e+1 2e+1];%2
    tau_I           = [1e+1 5e+1];%3
    p_pfu           = [8e+4 3e+5];%4
    p_rna           = [5e+6 8e+6];%5   
    tau_V           = [2e+1 5e+1];%6
    c_rna           = [5e-3 1e-2];%7
    omega_0         = [1e+0 5e+0];%8
    V0_pfu_stock    = [1e+5 3e+5];%9
    V0_rna_stock    = [1e+7 3e+7];%10
    V0_rna          = [1e+7 3e+7];%11
    n_L             = [1e+1 3e+1];%12
    n_V             = [1e+1 2e+1];%13

    pars = {beta;...
            tau_L;...
            tau_I;...
            p_pfu;...
            p_rna;...
            tau_V;...
            c_rna;...
            omega_0;...
            V0_pfu_stock;...
            V0_rna_stock;...
            V0_rna;...
            n_L;...
            n_V};
    
    sample=[];
    for ii=1:length(pars)
        a=log10(pars{ii}(1));
        b=log10(pars{ii}(2));
        step=(b-a)/n;
        interval_lower=a:step:b-step;
        interval_upper=a+step:step:b;
        sample=[sample;unifrnd(interval_lower,interval_upper,1,length(interval_lower))];
    end
 
    permparameters={};
    for ii=1:size(sample,1)
        iiparameter=sample(ii,:);
        iiparameter=iiparameter(randperm(length(iiparameter)));
        sample(ii,:) = iiparameter;
    end
    
    sample=sample';  
    save('paramsample.mat','sample');
    
end