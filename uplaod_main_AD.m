clear all
clc
close all
load data.mat
opts.alpha=0.01;
opts.r_e=0.50;
opts.r_mi=0.05;
opts.r_R=0.25;

opts.beta_e=0.25;
opts.beta_mi=0.25;
opts.beta_R=0.25;
neighbor=6;
[accuracy,spec,sens,PPV,NPV,AUC]=func_HGSCCA(train_epi,test_epi,train_miRNA,test_miRNA,train_RNA,test_RNA,train_label,test_label,opts,neighbor);


