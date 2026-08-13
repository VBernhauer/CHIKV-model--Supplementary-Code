function [] = linreg()

    clc;
    close all;

    %%% Vero cells IOL %%%
    data_IOL_RNA_Vero_Decay = readtable('../../Supplementary_Data.xlsx','Sheet','IOL RNA Vero Decay');

    time_Vero = [data_IOL_RNA_Vero_Decay.Var1;data_IOL_RNA_Vero_Decay.Var1;data_IOL_RNA_Vero_Decay.Var1];
    logV_Vero = log10([data_IOL_RNA_Vero_Decay.Var2;data_IOL_RNA_Vero_Decay.Var3;data_IOL_RNA_Vero_Decay.Var4]);

    mdl_Vero = fitlm(time_Vero,logV_Vero,'linear');
    disp(mdl_Vero);

end