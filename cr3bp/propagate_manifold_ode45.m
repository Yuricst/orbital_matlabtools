function [rr, vv, time] = ...
    propagate_manifold_ode45(mu,X0,Tmax,varargin)
% Function to propagate manifold with ode45()
% FORMAT: [rr, vv, time] = ...
%               propagate_manifold_ode45(mu,X0,Tmax,varargin)
% SIGN of Tmax should be NEGATIVE for stable manifold!
% varargin in order are:
%   stopaxis : 'x' or 'y' or 'z'
%   axisval : float to terminante manifold
%   nsteps
%   relTol
%   absTol
% Yuri Shimane, 2020/04/09
% ========================================================== %

% unpack varargin
if nargin == 3
    nsteps = 4000;
	relTol = 1e-7;
	absTol = 1e-8; 
elseif nargin == 5
    stopaxis = varargin{1};
    axisval = varargin{2};
    nsteps = 4000;
	relTol = 1e-7;
	absTol = 1e-8;
elseif nargin == 6
    stopaxis = varargin{1};
    axisval = varargin{2};
    nsteps = varargin{3};
	relTol = 1e-7;
	absTol = 1e-8;
elseif nargin ==7
    stopaxis = varargin{1};
    axisval = varargin{2};
    nsteps = varargin{3};
	relTol = varargin{4};
	absTol = 1e-8;
elseif nargin == 8
    stopaxis = varargin{1};
    axisval = varargin{2};
    nsteps = varargin{3};
	relTol = varargin{4};
	absTol = varargin{5};
end

% time array
time_ode45 = linspace(0,Tmax,nsteps)';
% format state input from cell to double
y0(1,1) = X0(1,1);  y0(2,1) = X0(2,1);  y0(3,1) = X0(3,1);
y0(4,1) = X0(4,1);  y0(5,1) = X0(5,1);  y0(6,1) = X0(6,1);
% call ODE45
if nargin == 3
    options = odeset('InitialStep',nsteps,...
            'RelTol',relTol,'AbsTol',absTol);
elseif nargin >= 5
    options = odeset('InitialStep',nsteps,'Events',@event_stopmanif,...
            'RelTol',relTol,'AbsTol',absTol);
end
[time,y] = ode45(@ode45_CR3BP,time_ode45,y0,options);  % fnc handle nested below
% extract position
rr = zeros(length(time),3);
rr(:,1) = y(:,1);
rr(:,2) = y(:,2)';
rr(:,3) = y(:,3)';
% extract velocity
vv = zeros(length(time),3);
vv(:,1) = y(:,4)';
vv(:,2) = y(:,5)';
vv(:,3) = y(:,6)';


% ================== Nested Functions ================== %
    % Event function for stopping criteria with 'on' Tf option
    function [value,isterminal,direction] = event_stopmanif(t, y)
        % when the y-value changes sign (=0)
        if strcmp(stopaxis,'x')
            value = y(1,1) - axisval;  % when value = 0
        elseif strcmp(stopaxis,'y')
            value = y(2,1) - axisval;  % when value = 0
        elseif strcmp(stopaxis,'z')
            value = y(3,1) - axisval;  % when value = 0
        end
        isterminal = 1;
        direction = 0;        
    end
    
    % CR3BP Equation of Motion
	function dydt = ode45_CR3BP(t,y)
        % radii in the CR3BP
        r1 = sqrt( (y(1,1)+mu)^2 + y(2,1)^2 + y(3,1)^2 );
        r2 = sqrt( (y(1,1)-1+mu)^2 + y(2,1)^2 + y(3,1)^2 );
        % velocities
        dydt(1,1) = y(4,1);
        dydt(2,1) = y(5,1);
        dydt(3,1) = y(6,1);
        % accelerations
        dydt(4,1) = 2*y(5,1) + y(1,1) - ((1-mu)/r1^3)*(mu+y(1,1)) + (mu/r2^3)*(1-mu-y(1,1));
        dydt(5,1) = -2*y(4,1) + y(2,1) - ((1-mu)/r1^3)*y(2,1) - (mu/r2^3)*y(2,1);
        dydt(6,1) = -((1-mu)/r1^3)*y(3,1) - (mu/r2^3)*y(3,1);
    end

end
