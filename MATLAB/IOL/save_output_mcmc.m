function [] = save_output_mcmc()
    %%% command:
    %%% save_output_mcmc()

    tmax        = 121;
    tt_decay    = 0:0.1:120;
    tt_infec    = 0:0.1:121;
    T0          = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fidm_decay_pfu  = fopen('output_mcmc_IOL_Vero_decay_pfu.txt','wt');
    fidm_decay_rna  = fopen('output_mcmc_IOL_Vero_decay_rna.txt','wt');
    fidm_pfu        = fopen('output_mcmc_IOL_Vero_pfu.txt','wt');
    fidm_rna        = fopen('output_mcmc_IOL_Vero_rna.txt','wt');

    T0_pfu  = 3e+5;
    MOI     = 0.1;
    V0_pfu  = MOI*T0_pfu;

    for rr = 1:5
        chains = load(strcat('./chains/chains_',num2str(rr),'.mat'));
        chains = chains.chains(:,:);
        for nn = 1:size(chains,1)
        
            if mod(nn,100) == 0
                disp([rr,nn]);
            end

            [logVdecay_out] = decay_Vero(10.^chains(nn,:));
            [logV_out]      = model_Vero(10.^chains(nn,:));

            fprintf(fidm_decay_pfu,'%.4f ',logVdecay_out{1});
            fprintf(fidm_decay_pfu,'\n');

            fprintf(fidm_decay_rna,'%.4f ',logVdecay_out{2});
            fprintf(fidm_decay_rna,'\n');

            fprintf(fidm_pfu,'%.4f ',logV_out{1});
            fprintf(fidm_pfu,'\n');

            fprintf(fidm_rna,'%.4f ',logV_out{2});
            fprintf(fidm_rna,'\n');
            
        end
    end
    fclose(fidm_decay_pfu);
    fclose(fidm_decay_rna);
    fclose(fidm_pfu);
    fclose(fidm_rna);
        

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
            V_pfu_ii = par(9)*(par(13)/par(6)*tt_decay).^(ii-1)/factorial(ii-1).*exp(-par(13)/par(6)*tt_decay);
            V_pfu = [V_pfu; V_pfu_ii];
        end
        logV_pfu = log10(sum(V_pfu,1));

        logV_rna = log10(par(10)*exp(-par(7)*tt_decay));

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
    function [logV] = model_Vero(par) 

        par(12) = round(par(12));
        par(13) = round(par(13));

        sol = ode23s(@ode,[0 tmax],[T0 zeros(1,par(12)) 0 V0_pfu zeros(1,par(13)-1) par(11)],[],par);
        
        PFU = [];        
        for kk = 3+par(12):2+par(12)+par(13)
            V_kk = deval(sol,tt_infec,kk);
            PFU = [PFU; V_kk];
        end
        logPFU = log10(sum(PFU,1));
        logRNA = log10(deval(sol,tt_infec,3+par(12)+par(13)));

        logV={logPFU;logRNA};

    end
         

end