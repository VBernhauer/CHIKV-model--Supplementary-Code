function [] = linreg()

    data_CARRIB_RNA_Decay = readtable('../../Supplementary_Data.xlsx','Sheet','CARRIB RNA Vero Decay');

    time = [data_CARRIB_RNA_Decay.Var1;data_CARRIB_RNA_Decay.Var1;data_CARRIB_RNA_Decay.Var1];
    logV = log10([data_CARRIB_RNA_Decay.Var2;data_CARRIB_RNA_Decay.Var3;data_CARRIB_RNA_Decay.Var4]);

    mdl = fitlm(time,logV,'linear');
    disp(mdl);

end