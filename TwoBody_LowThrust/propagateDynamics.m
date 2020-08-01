function [fval,c,ceq,Tfmax_update] = propagateDynamics(mu,u,Thrust,mdot,h0,hf,incf)
% Functon propagates dynamics with ode45

% unpack control vector
nsteps = (length(u)-1)/3;
tof = u(1,1);   % indx 1
umag = u(2:1+nsteps,1);  % indx 2 ~ 101
ualfa = u(2+nsteps:1+2*nsteps,1);  % indx 102 ~ 201
ubeta = u(2+2*nsteps:1+3*nsteps,1);  % indx 202 ~ 301 (end)

% set time-span
tspan = linspace(0,tof,nsteps);

% construct initial state: [x, y, z, vx, vy, vz, mass]
X0 = [];

% call ode45
odeopt = odeset('Events',@eventReachOrbit);
[time,dynmat] = ode45(@twobody,tspan,X0,odeopt);

% decompose final state
Xf = dynmat(end,:);
% final altitude
rp = sqrt(Xf(1,1)^2 + Xf(1,2)^2 + Xf(1,3)^2);
vp = sqrt(Xf(1,4)^2 + Xf(1,5)^2 + Xf(1,6)^2);
hf - cross(rp, vp);
incp = acosd(hf(3) / norm(hf));

% compute objective function
fval = -dynmat(end,7);  % maximize final mass

% compute non-linear inequality constraints
c = [];

% compute non-linear equality constraints
ceq(1) = rp - hf;  % altitude constraint (circularize)
ceq(2) = vp - sqrt(1/hf);  % velocity constraint (circularize)
ceq(3) = incp - incf;      % inclination

% ... Nested functions ... %
    % Equation of motion
    function Xdot = twobody(t,X)
        % unpack state
        x = X(1,1); 
        y = X(2,1); 
        z = X(3,1);
        vx = X(4,1);
        vy = X(5,1);
        vz = X(6,1);
        m = X(7,1);
        % compute radius
        r = sqrt(x^2 + y^2 + z^2);
        % interpolate control parameter at time t
        umag_t = interp1(tspan, umag, t);
        ualfa_t = interp1(tspan, ualfa, t);
        ubeta_t = interp1(tspan, ubeta, t);
        % compute thrust acceleration
        ax_Thrust = umag_t * cosd() * sind();
        ay_Thrust = umag_t;
        az_Thrust = umag_t;
        % position derivative (= velocity)
        Xdot(1,1) = vx;
        Xdot(2,1) = vy;
        Xdot(3,1) = vz;
        % velocity derivative (= acceleration)
        Xdot(4,1) = -(mu/r^3)*x + ax_Thrust;
        Xdot(5,1) = -(mu/r^3)*y + ay_Thrust;
        Xdot(6,1) = -(mu/r^3)*z + az_Thrust;
        % mass derivative
        Xdot(7,1) = mdot;
    end

%     % Event function (ode45 stop condition)
%     function [altitude,isterminal,direction] = eventReachOrbit(t,X)
%         % unpack state
%         x = X(1,1);
%         y = X(2,1);
%         z = X(3,1);
%         % compute norm
%         r = sqrt(x^2 + y^2 + z^2);
%         altitude = hf - r; % The value that we want to be zero
%         isterminal = 1;    % Halt integration 
%         direction = 0;     % The zero can be approached from either direction
%     end

% ======================== %
end