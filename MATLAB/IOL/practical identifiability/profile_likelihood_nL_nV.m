function [results] = profile_likelihood_nL_nV(nstart,nend,n_L_span)
    
    clc;
    close all;
    rng('shuffle');

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
    T0      = 1;
    T0_pfu  = 3e+5;
    MOI     = 0.1;
    V0_pfu  = MOI*T0_pfu;

    tmax        = 121;
    time_decay  = [0,2,4,8,16,24,48,72,96,120];
    time_infec  = [1.5,3,5,9,17,25,49,73,97,121];

    MLpars_table    = readmatrix('../maxLikValues.txt');
    MLpars          = transpose(MLpars_table(2:end,2));
    MLvalue         = transpose(MLpars_table(1,2));
    MLvalues        = [MLvalue,MLpars];
    save('./MLvalues.mat','MLvalues');

    n_V_span  = 12:1:24;

    PopSz = 5;   

    for iround=nstart:nend

        initPopul = [];
        for ll = 1:11
            parRnd = MLpars(ll)*0.975 + (MLpars(ll)*1.025 - MLpars(ll)*0.975).*rand(PopSz,1);
            initPopul = [initPopul, parRnd];
        end          
        funval = [];  
        for j=1:length(n_V_span)
            idx_est = [1:11];
            idx_fix = [12,13];
            options = optimoptions(@fmincon,'Display','off','Algorithm','interior-point');
            fval_opt = 9999.99;
            for mm = 1:size(initPopul,1)
                [x,fval,exitflag,output] = fmincon(@(params)loglike(params,[n_L_span(iround),n_V_span(j)],idx_est,idx_fix), initPopul(mm,1:11), [], [], [], [], 0*ones(length(idx_est),1), inf*ones(length(idx_est),1), [], options); 
                if fval < fval_opt
                    fval_opt = fval;                        
                end
            end
            funval = [funval; n_L_span(iround), n_V_span(j), fval_opt];            
        end  
        results{iround} = funval;
        save(strcat('./PL_nL_nV/nL_',num2str(n_L_span(iround)),'_nV.mat'),'funval'); 

    end

    disp('Done.')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% helper functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% log-likelihood formula %%%
    function [value] = lognormpdf(data,model,stdev)
        value = (-1)*(-0.5*((data-model)./stdev).^2  - log(sqrt(2*pi).*stdev));
    end

    %%% log-likelihood %%%
    function [value] = loglike(pars_est,pars_fix,idx_est,idx_fix) 
        logV_decay  = decay(pars_est,pars_fix,idx_est,idx_fix);
        logV_inf    = infection(pars_est,pars_fix,idx_est,idx_fix);
        logV        = {logV_decay;logV_inf};
        
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

    %%% viral decay %%%
    function [output] = decay(pars_est,pars_fix,idx_est,idx_fix) 
        [beta,...
        tau_L,...
        tau_I,...
        p_pfu,...
        p_rna,...
        tau_V,...
        c_rna,...
        omega_0,...
        V0_pfu_stock,...
        V0_rna_stock,...
        V0_rna,...
        n_L,...
        n_V]=findpars(pars_est,pars_fix,idx_est,idx_fix);

        n_V = round(n_V);

        V_pfu = [];
        for ii = 1:n_V
            V_pfu_ii = V0_pfu_stock*exp(-n_V/tau_V*time_decay).*((n_V/tau_V*time_decay).^(ii-1))/factorial(ii-1);
            V_pfu = [V_pfu; V_pfu_ii];
        end
        logV_pfu = log10(sum(V_pfu,1));
        logV_rna = log10(V0_rna_stock*exp(-c_rna*time_decay));

        output = {logV_pfu;logV_rna};
    end

    %%% output of the infection model %%%
    function [output] = infection(pars_est,pars_fix,idx_est,idx_fix)
        [beta,...
        tau_L,...
        tau_I,...
        p_pfu,...
        p_rna,...
        tau_V,...
        c_rna,...
        omega_0,...
        V0_pfu_stock,...
        V0_rna_stock,...
        V0_rna,...
        n_L,...
        n_V]=findpars(pars_est,pars_fix,idx_est,idx_fix);        

        n_L = round(n_L);
        n_V = round(n_V);

        sol = ode23s(@ode,[0 tmax],[T0 zeros(1,n_L) 0 V0_pfu zeros(1,n_V-1) V0_rna], [], pars_est, pars_fix, idx_est, idx_fix);
        
        PFU = [];        
        for kk = 3+n_L:2+n_L+n_V
            V_kk = deval(sol,time_infec,kk);
            PFU = [PFU; V_kk];
        end
        logPFU = log10(sum(PFU,1));
        logRNA = log10(deval(sol,time_infec,3+n_L+n_V));

        output={logPFU;logRNA};
    end

    %%% ODE system %%%
    function [xdot] = ode(t,x,pars_est,pars_fix,idx_est,idx_fix)     
       [beta,...
        tau_L,...
        tau_I,...
        p_pfu,...
        p_rna,...
        tau_V,...
        c_rna,...
        omega_0,...
        V0_pfu_stock,...
        V0_rna_stock,...
        V0_rna,...
        n_L,...
        n_V]=findpars(pars_est,pars_fix,idx_est,idx_fix);

        n_L = round(n_L);
        n_V = round(n_V);

        xdot = zeros(3+n_L+n_V,1);
        
        %%% washing %%%
        omega = washing(t,omega_0);

        %%% model %%%
        xdot(1) = -beta*(sum(x(3+n_L:2+n_L+n_V)))*x(1);    
        xdot(2) = beta*(sum(x(3+n_L:2+n_L+n_V)))*x(1) - n_L/tau_L*x(2);
        for ii = 3:1+n_L
            xdot(ii) = n_L/tau_L*(x(ii-1) - x(ii));
        end
        xdot(2+n_L) = n_L/tau_L*x(1+n_L) - 1/tau_I*x(2+n_L);
        xdot(3+n_L) = p_pfu*x(2+n_L) - omega*x(3+n_L) - n_V/tau_V*x(3+n_L);
        for kk = 4+n_L:2+n_L+n_V
           xdot(kk) = n_V/tau_V*(x(kk-1) - x(kk)) - omega*x(kk);
        end      
        xdot(3+n_L+n_V) = p_rna*x(2+n_L) - omega*x(3+n_L+n_V) - c_rna*x(3+n_L+n_V);
    end

    %%% washing %%%
    function w = washing(t,w0)

        % w0 ... strength of washing
        % wd ... standard deviation of the length of washing (hours)
        % wt ... time of washing implementation (hours)
        wd = 0.05;
        wt = 1; 
        w = w0*1/sqrt(2*pi*wd^2)*exp(-(t-wt)^2/(2*wd^2));

    end

    %%% fixed and estimated parameters %%%
    function [beta,...
              tau_L,...
              tau_I,...
              p_pfu,...
              p_rna,...
              tau_V,...
              c_rna,...
              omega_0,...
              V0_pfu_stock,...
              V0_rna_stock,...
              V0_rna,...
              n_L,...
              n_V]=findpars(pars_est,pars_fix,idx_est,idx_fix)
    
        if length(find(idx_fix==1))==0
            I=find(idx_est==1);
            beta = pars_est(I);
        else
            I=find(idx_fix==1);
            beta = pars_fix(I);
        end
        
        if length(find(idx_fix==2))==0
            I=find(idx_est==2);
            tau_L = pars_est(I);
        else
            I=find(idx_fix==2);
            tau_L = pars_fix(I);
        end
        
        if length(find(idx_fix==3))==0
            I=find(idx_est==3);
            tau_I = pars_est(I);
        else
            I=find(idx_fix==3);
            tau_I = pars_fix(I);
        end
        
        if length(find(idx_fix==4))==0
            I=find(idx_est==4);
            p_pfu = pars_est(I);
        else
            I=find(idx_fix==4);
            p_pfu = pars_fix(I);
        end     

        if length(find(idx_fix==5))==0
            I=find(idx_est==5);
            p_rna = pars_est(I);
        else
            I=find(idx_fix==5);
            p_rna = pars_fix(I);
        end   

        if length(find(idx_fix==6))==0
            I=find(idx_est==6);
            tau_V = pars_est(I);
        else
            I=find(idx_fix==6);
            tau_V = pars_fix(I);
        end   

        if length(find(idx_fix==7))==0
            I=find(idx_est==7);
            c_rna = pars_est(I);
        else
            I=find(idx_fix==7);
            c_rna = pars_fix(I);
        end   

        if length(find(idx_fix==8))==0
            I=find(idx_est==8);
            omega_0 = pars_est(I);
        else
            I=find(idx_fix==8);
            omega_0 = pars_fix(I);
        end   

        if length(find(idx_fix==9))==0
            I=find(idx_est==9);
            V0_pfu_stock = pars_est(I);
        else
            I=find(idx_fix==9);
            V0_pfu_stock = pars_fix(I);
        end   

        if length(find(idx_fix==10))==0
            I=find(idx_est==10);
            V0_rna_stock = pars_est(I);
        else
            I=find(idx_fix==10);
            V0_rna_stock = pars_fix(I);
        end   

        if length(find(idx_fix==11))==0
            I=find(idx_est==11);
            V0_rna = pars_est(I);
        else
            I=find(idx_fix==11);
            V0_rna = pars_fix(I);
        end   

        if length(find(idx_fix==12))==0
            I=find(idx_est==12);
            n_L = pars_est(I);
        else
            I=find(idx_fix==12);
            n_L = pars_fix(I);
        end   

        if length(find(idx_fix==13))==0
            I=find(idx_est==13);
            n_V = pars_est(I);
        else
            I=find(idx_fix==13);
            n_V = pars_fix(I);
        end   
       
    end



end