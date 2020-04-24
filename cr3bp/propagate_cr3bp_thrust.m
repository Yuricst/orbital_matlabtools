function [time,dynmat] = propagate_cr3bp_thrust(mu,u,mass0,X0,Isp,Fmax,nsteps)
% Function propagates dynamics of orbit in CR3BP and 
% compute objective function value and nonlinear conditions
% Yuri Shimane, 2020/04/17
% ============================================================== %

% length of time-steps
timestep = (length(u)-1)/3;

% unpack state to get final time
ux = u(1:timestep, 1);
uy = u(timestep+1:2*timestep , 1);
uz = u(2*timestep+1:3*timestep , 1);
Tf = u(end);
% time array
time_cntrl = linspace(0,Tf,timestep);

% include mass to X0
X0mass = X0;
X0mass(7, 1) = mass0;

% setup option for ode45 to integrate dynamics
relTol = 1e-8;
absTol = 1e-10;
opts = odeset('InitialStep',nsteps,'RelTol',relTol,'AbsTol',absTol);

% integrade dynamics
[time, dynmat] = ode45(@(t,X) rhs_cr3bp(t, X, time_cntrl, ux, uy, uz),[0 Tf],X0mass,opts);

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
        % compute thrust magnitude
        T = Fmax/mass;
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
        Xdot(4,1) = 2*vy + x - ((1-mu)/r1^3)*(mu+x) + (mu/r2^3)*(1-mu-x) + ux_val*T;
        Xdot(5,1) = -2*vx + y - ((1-mu)/r1^3)*y - (mu/r2^3)*y + uy_val*T;
        Xdot(6,1) = -((1-mu)/r1^3)*z - (mu/r2^3)*z + uz_val*T;
        % --- MASS DERIVATIVE --- %
        Xdot(7,1) = -Fmax/(Isp*9.81);
    end

end