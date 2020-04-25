% testing lambert 
% earth-to-mars 2d transfer

% house keeping
clear; close all; clc;
% test case (Earth-Mars transfer)
mu = 132712000000;
tof_days = 210;         % [days]
tof = tof_days*24*60*60; % [sec]
AU = 149597871;        % 1 AU in km
phi = 40;              % approximate phasing angle for hohmann transfer [deg]
% initial state
r0 = [1*AU*cosd(phi),1*AU*sind(phi),0];
% final state
rf = [-1.52*AU,0,0];

% solve lambert
[v0,vf] = lambert_curtis(r0,rf,tof,mu,'prograde');

