function [] = save_output_mcmc()

    close all
    set(0,'DefaultFigureVisible','on');

    tmax        = 73;
    tt_decay    = 0:0.1:72;
    tt_infec    = 0:0.1:73;
    T0          = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fidm_decay_pfu  = fopen('output_mcmc_CARRIB_decay_pfu.txt','wt');
    fidm_decay_rna  = fopen('output_mcmc_CARRIB_decay_rna.txt','wt');
    fidm_pfu_001    = fopen('output_mcmc_CARRIB_pfu_001.txt','wt');
    fidm_rna_001    = fopen('output_mcmc_CARRIB_rna_001.txt','wt');
    fidm_pfu_1      = fopen('output_mcmc_CARRIB_pfu_1.txt','wt');
    fidm_rna_1      = fopen('output_mcmc_CARRIB_rna_1.txt','wt');

    for rr = 1:5
        chains = load(strcat('./chains/chains_',num2str(rr),'.mat'));
        chains = chains.chains(:,:);
        for nn = 1:size(chains,1)
            
            if mod(nn,100) == 0
                disp([rr,nn]);
            end

            [decay_out] = decay(10.^chains(nn,:));
            [logV_out] = model(10.^chains(nn,:));

            fprintf(fidm_decay_pfu,'%.4f ',decay_out{1});
            fprintf(fidm_decay_pfu,'\n');

            fprintf(fidm_decay_rna,'%.4f ',decay_out{2});
            fprintf(fidm_decay_rna,'\n');

            fprintf(fidm_pfu_001,'%.4f ',logV_out{1}{1});
            fprintf(fidm_pfu_001,'\n');

            fprintf(fidm_rna_001,'%.4f ',logV_out{1}{2});
            fprintf(fidm_rna_001,'\n');

            fprintf(fidm_pfu_1,'%.4f ',logV_out{2}{1});
            fprintf(fidm_pfu_1,'\n');

            fprintf(fidm_rna_1,'%.4f ',logV_out{2}{2});
            fprintf(fidm_rna_1,'\n');
        end
    end

    fclose(fidm_decay_pfu);
    fclose(fidm_decay_rna);
    fclose(fidm_pfu_001);
    fclose(fidm_rna_001);
    fclose(fidm_pfu_1);
    fclose(fidm_rna_1);


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

        logV_pfu = log10(par(8)*exp(-par(6)*tt_decay));
        logV_rna = log10(par(9)*ones(1,length(tt_decay)));

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

        %%% low MOI        
        sol_low = ode23s(@ode,[0 tmax],[T0 zeros(1,par(10)) 0 0.01*par(8) 0.01*par(9)],[],par);
        logPFU_low = log10(deval(sol_low,tt_infec,3+par(10)));
        logRNA_low = log10(deval(sol_low,tt_infec,4+par(10)));

        %%% high MOI
        sol_high = ode23s(@ode,[0 tmax],[T0 zeros(1,par(10)) 0 par(8) par(9)],[],par);
        logPFU_high = log10(deval(sol_high,tt_infec,3+par(10)));
        logRNA_high = log10(deval(sol_high,tt_infec,4+par(10)));

        logV={{logPFU_low;logRNA_low},{logPFU_high;logRNA_high}};

    end
         

end