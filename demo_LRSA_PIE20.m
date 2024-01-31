close all; clear all;clc 
addpath(genpath('./'));
%% Data Preparation 
load('PIE_20.mat','X','gnd');
nGroup=20;
%% parameters setting
    lambda = 1;  %here lambda is used to deal with the loss term
    beta = 600;      %here beta is used to balance the low-rank and sparse terms
    
%% Main algorithm
    W=LRSA(X,lambda,beta);
    
%% Different evaluation metrics  
    INI=[10,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,390];   % INI is used to set the cluster center
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
    
   save('PIE20LRSA.mat','measure')  