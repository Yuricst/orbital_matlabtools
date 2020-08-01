function dynmat = dynamics_orbit_raising(time,u)
% function propagates the dynamics of orbit raising (integrator)
% INPUT
%   time : nsteps x 1 time array
%   u : 2*nsteps x 1 control input array

nsteps = length(time);
dt = time(1,2) - time(1,1);
% handle for rhs expressions of dynamical system
rhs = @rhs_orbit_raising;

% initialize array to store
r      = zeros(nsteps,1);
theta  = zeros(nsteps,1);
vr     = zeros(nsteps,1);
vtheta = zeros(nsteps,1);
% initial conditions
r(1,1)      = 1;
theta(1,1)  = 0;
vr(1,1)     = 0;
vtheta(1,1) = 1;
% initialize waitbar
% f = waitbar(0,'Propagating dynamics....','Name','Progress');
% integrate rhs equation
for i = 1:nsteps-1
    % update waitbar
%     waitbar(i/(nsteps-1),f,['Propagating dynamics (',num2str(round(i/(nsteps-1)*100)),'%)']);
    % compute rhs
    [rdot,thetadot,vrdot,vthetadot] = rhs(r(i,1),theta(i,1),vr(i,1),...
        vtheta(i,1),u(i,1),u(nsteps+i,1),time(i));
    % update with time-step (EULER method)
    r(i+1,1)      = r(i,1)      + dt*rdot;
    theta(i+1,1)  = theta(i,1)  + dt*thetadot;
    vr(i+1,1)     = vr(i,1)     + dt*vrdot;
    vtheta(i+1,1) = vtheta(i,1) + dt*vthetadot;
end
% close waitbar
% close(f);

dynmat = [r, theta, vr, vtheta];

end

