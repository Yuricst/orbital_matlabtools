function [fval,gradf,c,ceq,gradc,gradceq] = rundynamics(u,nsteps,tf,mcxop)
% function computes objective function, non-linear constraints
% =================================================== %
% Yuri Shimane, 2020/02/20
% INPUT
%   u : control input
%   nsteps : number of steps
%   tf : final time
%   mcxop : string whether to compute gradient with MCX vectors
% =================================================== %

% create time array
time = linspace(0,tf,nsteps);
% integrate dynamics
if mcxop == false
    dynmat = dynamics_orbit_raising(time,u);
    gradf = [];
elseif mcxop == true
    % gradient is also returned from MCX integration
    [dynmat,gradf] = dynamics_orbit_raising_MCX(time,u);
end

% objective function value: (negative of) altitude
fval = -dynmat(nsteps,1);
% plot_evol_orbit(time,dynmat);

% non-linear inequality constraints
c = zeros(1,nsteps);
for i = 1:nsteps
    % relaxed path constraint
    c(1,i) = u(i,1)^2 + u(nsteps+i,1)^2 - 1;
end

% non-linear equality constraints
% ceq(1,1) and ceq(2,1) : initial and terminal boundary conditions
ceq(1,1) = dynmat(1,1) - 1;  % r(0) = 1
ceq(1,2) = dynmat(1,2);      % theta(0) = 0
ceq(1,3) = dynmat(1,3);      % vr(0) = 0
ceq(1,4) = dynmat(1,4) - 1;  % vtheta(0) = 1
ceq(1,5) = dynmat(nsteps,3); % vr(tf) = 0
ceq(1,6) = sqrt(1/dynmat(nsteps,1)) - dynmat(nsteps,4); % sqrt(1/r(tf)) - vtheta(tf) = 0
% % ceq(1,7:end) : path constraint
% for i = 1:nsteps
%     ceq(1,i+6) = u(i,1)^2 + u(nsteps+i,1)^2 - 1;
% end

% sensitivity to non-linear conditions
gradc = 0;
gradceq = 0;

end
