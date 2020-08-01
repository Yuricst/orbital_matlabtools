function [rdot,thetadot,vrdot,vthetadot] = rhs_orbit_raising(r,theta,vr,vtheta,u1,u2,t)
% define right-hand side of dynamic system for orbit raising problem
% Yuri Shimane, 2020/02/07
% =========================================================== %
% INPUT:
%   states : states in order r, theta, vr, vtheta, u1, u2
%   t : time after starting to integrate
% =========================================================== %

% compute acceleration
a = 0.1405/(1-0.0749*t);

% dynamic system equations
rdot      = vr;
thetadot  = vtheta/r;
vrdot     = vtheta^2/r - 1/r^2 + a*u1;
vthetadot = -vtheta*vr/r + a*u2;

end

