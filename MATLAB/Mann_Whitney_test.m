function [] = Mann_Whitney_test()

    %%% command:
    %%% Mann_Whitney_test()
   
    clc;
    close all
    set(0,'DefaultFigureVisible','on');

    carrib_chains = [];
    for jj = 1:5
        jjchains = load(strcat('./CARRIB/chains/chains_',num2str(jj),'.mat'));
        jjchains = jjchains.chains(:,:);
        for kk = 1:size(jjchains,1)
            carrib_chains = [carrib_chains;jjchains(kk,:)];
        end
    end

    iol_chains = [];
    for jj = 1:5
        jjchains = load(strcat('./IOL/chains/chains_',num2str(jj),'.mat'));
        jjchains = jjchains.chains(:,:);
        for kk = 1:size(jjchains,1)
            iol_chains = [iol_chains;jjchains(kk,:)];
        end
    end

    %%% beta
    [p,h,stats] = ranksum(carrib_chains(:,1), iol_chains(:,1))

    %%% tau_L
    [p,h,stats] = ranksum(carrib_chains(:,2), iol_chains(:,2))

    %%% tau_I
    [p,h,stats] = ranksum(carrib_chains(:,3), iol_chains(:,3))

    %%% p_pfu
    [p,h,stats] = ranksum(carrib_chains(:,4), iol_chains(:,4))

    %%% p_rna
    [p,h,stats] = ranksum(carrib_chains(:,5), iol_chains(:,5))

    %%% omega_0
    [p,h,stats] = ranksum(carrib_chains(:,7), iol_chains(:,8))

    %%% n_L
    [p,h,stats] = ranksum(round(10.^carrib_chains(:,10)), round(10.^iol_chains(:,13)))

end