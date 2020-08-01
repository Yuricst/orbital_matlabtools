%% main - low-thrust transfer from LEO to GEO

% house keep
clear; close all; clc;

%% setup problem, initial guess, optim options
% initial and target orbit
RE = 6378;  % [km]
h0 = RE + 200;  % [km]
inc0 = 0; % [deg]
hf = RE + 35786; % [km]
incf = 0;  % [deg]
mu = 398600;

% thruster
Isp = 1500;
mdot = -40;
Thrust = 0.015;  % [N]

% initial guess for control vector: Tf, ux/uy/uz or u, alpha, beta
u0_struct.tof = 30*24*60*60; % guess on Tf [sec]
u0_struct.u = zeros(100,1);
u0_struct.alfa = zeros(100,1);
u0_struct.beta = zeros(100,1);

% run function calling optimizer
[x,fval,exitflag,output] = runoptim(mu,Thrust,u0_struct,h0,hf,incf,mdot);

%   in the optimizer...
%       - propagate dynamics until some final time Tf
%       - constraint: circularize with correct inclination







