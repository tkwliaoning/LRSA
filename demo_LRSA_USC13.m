close all; clear all;clc 
addpath(genpath('./'));
%% Data Preparation 
load('USC13_13¡Á80_33112_lbp_riu2.mat','X','gt');
nGroup=13;
gnd=gt;
%% parameters setting
    lambda = 6200;  %here lambda is used to deal with the loss term
    beta = 10;      %here beta is used to balance the low-rank and sparse terms
    
%% Main algorithm
    W=LRSA(X,lambda,beta);
    
%% Different evaluation metrics  
    INI=[1:80:961];   % INI is used to set the cluster center
    label = fixSpectralClustering(W,nGroup,INI);
    result = ClusteringMeasure(gnd, label);
    measure(1,1)=result(1);
    measure(2,1)=result(2);
    measure(3,1)=result(3);
    measure(4,1)=RandIndex(gnd, label);
    [fscore,p,r] = compute_f(gnd, label);
    measure(5,1)=fscore;
    measure(6,1)=p;
    measure(7,1)=r;
    measure(8,1)=lambda;
    measure(9,1)=beta;
    
   save('USC13LRSA.mat','measure')  