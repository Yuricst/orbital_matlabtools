function [fval,c,ceq] = rundynamics_raiseOrb(mu,u,mass0,X0,Xf,Isp,Fmax,nsteps,Lstar,Tstar)
% Function propagates dynamics of orbit in CR3BP and 
% compute objective function value and nonlinear conditions
% Yuri Shimane, 2020/04/17
% ============================================================== %

% unpack state to get final time
ux = u(1:nsteps, 1);
uy = u(nsteps+1:2*nsteps , 1);
uz = u(2*nsteps+1:3*nsteps , 1);
Tf = u(end);
% time array
time_cntrl = linspace(0,Tf,nsteps);

% include mass to X0
X0mass = X0;
X0mass(7, 1) = mass0;

% setup option for ode45 to integrate dynamics
relTol = 1e-8;
absTol = 1e-10;
opts = odeset('InitialStep',nsteps,'RelTol',relTol,'AbsTol',absTol);

% integrade dynamics
[time, dynmat] = ode45(@(t,X) rhs_cr3bp(t, X, time_cntrl, ux, uy, uz),time_cntrl,X0mass,opts);
% size of dynmat
[n,~] = size(dynmat);

% compute objective function value from dynmat
fval = dynmat(n,7);  % minimize mass

% non-linear inequality constraints (path control)
c = zeros(1,nsteps);
for i = 1:nsteps
    c(1,i) = ux(i,1)^2 + uy(i,1)^2 + uz(i,1)^2 - 1;
end

% non-linear equality constraints (match final state)
ceq = zeros(1,6);
for j = 1:6
    ceq(1,j) = dynmat(n,j) - Xf(j,1);
end

% ... nested functions below ...   
    % === RHS of dynamics === %
    function Xdot = rhs_cr3bp(t, X, time_cntrl, ux, uy, uz)
        % Function computes rhs of dynamics including mass as 7th state 
        % ====================================================== %
        % decompose state and STM
        x = X(1,1);  y = X(2,1);  z = X(3,1);
        vx = X(4,1); vy = X(5,1); vz = X(6,1);
        mass = X(end);
        
        % interpolate the control set (time,ux) at time t
        ux_val = interp1(time_cntrl,ux,t);
        % interpolate the control set (time,uy) at time t
        uy_val = interp1(time_cntrl,uy,t);
        % interpolate the control set (time,uz) at time t
        uz_val = interp1(time_cntrl,uz,t);
        % compute acceleration magnitude
        T = Fmax/mass;  % [m/sec^2]
        Tnondim = T * (Tstar^2/Lstar);
        % calculate radii
        r1 = sqrt( (x+mu)^2 + y^2 + z^2 );
        r2 = sqrt( (x-1+mu)^2 + y^2 + z^2 );
        % --- STATE DERIVATIVE --- %
        % initialize double
        Xdot = zeros(7,1);
        % position-state derivative
        Xdot(1,1) = vx;
        Xdot(2,1) = vy;
        Xdot(3,1) = vz;
        % velocity-state derivative
        Xdot(4,1) = 2*vy + x - ((1-mu)/r1^3)*(mu+x) + (mu/r2^3)*(1-mu-x) + ux_val*Tnondim;
        Xdot(5,1) = -2*vx + y - ((1-mu)/r1^3)*y - (mu/r2^3)*y + uy_val*Tnondim;
        Xdot(6,1) = -((1-mu)/r1^3)*z - (mu/r2^3)*z + uz_val*Tnondim;
        % --- MASS DERIVATIVE --- %
        Xdot(7,1) = -Fmax/(Isp*9.81);  % [kg/sec]
    end

end