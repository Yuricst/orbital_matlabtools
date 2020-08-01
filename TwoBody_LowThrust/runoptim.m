function [u,fval,exitflag,output] = runoptim(mu,Thrust,u0_struct,h0,hf,incf,mdot)
% Function calls fminconfor low-thrust orbit raising problem
% Yuri Shimane, 2020/05/07

% declare global variables for preventing unecessary dynamics computation
xLast = [];
fvalLast = [];
cLast = [];
ceqLast = [];

% non-dimensionalization parameter
r_star = h0;  % [km]
v_star = sqrt(mu/h0);  % [km/sec]
t_star = r_star / v_star; % [sec]
F_star = 1;  % [N]

% non-dimensionalize input
h0_nd = h0/r_star;
hf_nd = hf/r_star;
Thrust_nd = Thrust;

% unpack initial guess structure u0struct
u0_nd.tof = u0_struct.tof / t_star;  % (non-dim time)
u0_nd.u = u0_struct.u / F_star;      % (non-dim force)
u0 = vertcat(u0_nd.tof, u0_nd.u, u0_struct.alfa, u0_struct.beta);

% bounds on u0_nd
nsteps = length(u0_struct.u);
lb = [0; -ones(nsteps,1); -ones(nsteps,1); -ones(nsteps,1)];
ub = [u0_nd.tof*10; -ones(nsteps,1); -ones(nsteps,1); -ones(nsteps,1)];

% set options
opts = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',1.5e5,... 
        'MaxIterations',4e3,'ConstraintTolerance',1.0000e-03,...
        'FiniteDifferenceType','forward','FiniteDifferenceStepSize',1e-8);

% run fmincon
[u,fval,exitflag,output] =  fmincon(@objfun,u0,[],[],[],[],lb,ub,@nlcon,opts);

% propagate final solution
%[time,dynmat] = ode45(@twobody,[0 u(1,1)],X0);

% ... nested functions ... %
    
    % Objective function
    function fval = objfun(u)
        if ~isequal(u,xLast)
            % run dynamics
            [fvalLast, cLast, ceqLast] = propagateDynamics(mu,u,Thrust_nd,mdot,h0_nd,hf_nd,incf);
            xLast = u;
        end
        fval = fvalLast;
    end

    % Non-linear constraints
    function [c,ceq] = nlcon(u)
        if ~isequal(u,xLast)
            % run dynamics
            %[fvalLast, ...] = function()
            xLast = u;
        end
        c = cLast;
        ceq = ceqLast;
    end

    % Equation of motion for propagating final result
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
        % compute thrust acceleration
        anorm = Thrust/m;
        % position derivative (= velocity)
        Xdot(1,1) = vx;
        Xdot(2,1) = vy;
        Xdot(3,1) = vz;
        % velocity derivative (= acceleration)
        Xdot(4,1) = -(mu/r^3)*x + anorm;
        Xdot(5,1) = -(mu/r^3)*y + anorm;
        Xdot(6,1) = -(mu/r^3)*z + anorm;
        % mass derivative
        Xdot(7,1) = mdot;
    end
% ------------------------ % 
end
