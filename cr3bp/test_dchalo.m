%% Testing for Halo construction with DC

% house keep 
clear; close all; clc;

% analytical Richardson 3rd order solution 
X0_input = [0.82568901, 0.        , 0.07478352, 0.        , 0.19170699,       0.        ]';
period_richard = 2.7757924329382075;
Thalf = period_richard/2;

% define cr3bp system
mu = 0.012409319;
lstar = 384400;
residualtol = 1e-8;

[X0_halo,Tfull,varargout] = ...
    collinearHalo_differentialCorrection(mu,X0_input,Thalf,residualtol);




