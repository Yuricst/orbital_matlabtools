function [rr, vv, time] = ...
    propagate_manifold_poincare(mu,X0,Tmax,varargin)
% Function to propagate manifold with ode45() until prescribed poincare
% section
% FORMAT: [rr, vv, time] = ...
%               propagate_manifold_ode45(mu,X0,Tmax,varargin)
% SIGN of Tmax should be NEGATIVE for stable manifold!
% varargin options are the four Poincare sections
%   U1 = {(x,y,z); x<0, y=0}
%   U2 = {(x,y,z); x=1-mu, y<0}
%   U3 = {(x,y,z); x=1-mu, y>0}
%   U4 = {(x,y,z); x<-1, y=0}
% Yuri Shimane, 2020/04/25
% ========================================================== %

% unpack varargin
numarg = length(varargin);
% expect up to two poincare sections
if numarg == 1
    PS1 = varargin{1};
elseif numarg == 2
    PS1 = varargin{1};
    PS2 = varargin{2};
end

% setup for ODE45
nsteps = 4000;
relTol = 1e-3;
absTol = 1e-6; 
% time array
time_ode45 = linspace(0,Tmax,nsteps)';

% set ODE45 option
if numarg == 0
    options = odeset('InitialStep',nsteps,...
            'RelTol',relTol,'AbsTol',absTol);
elseif numarg >= 1
    options = odeset('InitialStep',nsteps,'Events',@eventPoincare,...
            'RelTol',relTol,'AbsTol',absTol);
end
% call ODE45
[time,dynmat] = ode45(@ode45_CR3BP,time_ode45,X0,options);  % fnc handle nested below
% extract position
rr = zeros(length(time),3);
rr(:,1) = dynmat(:,1);
rr(:,2) = dynmat(:,2);
rr(:,3) = dynmat(:,3);
% extract velocity
vv = zeros(length(time),3);
vv(:,1) = dynmat(:,4);
vv(:,2) = dynmat(:,5);
vv(:,3) = dynmat(:,6);


% ================== Nested Functions ================== %
    % Event function for stopping criteria with 'on' Tf option
    function [value,isterminal,direction] = eventPoincare(t, y)
        % check number of sections (expect up to 2)
        if strcmp(PS1,'U1')
            val1 = y(1,1)<0 && y(2,1)>=0; %[y(1,1)+1e-8, y(2,1)];
            % direction (depends on stable or unstable manifold)
            if Tmax < 0 % stable manifold
                dir1 = 1;
            else        % unstable manifold
                dir1 = -1;
            end
        elseif strcmp(PS1,'U2')
            val1 = y(1,1)<=1-mu && y(2,1)<0;
            % direction (depends on stable or unstable manifold)
            if Tmax < 0 % stable manifold
                dir1 = -1;
            else        % unstable manifold
                dir1 = 1;
            end
        elseif strcmp(PS1,'U3')
            val1 = y(1,1)>=1-mu && y(2,1)>0;
            % direction (depends on stable or unstable manifold)
            if Tmax < 0 % stable manifold
                dir1 = 1;
            else        % unstable manifold
                dir1 = -1;
            end
        elseif strcmp(PS1,'U4')
            val1 = y(1,1)<-1 && y(2,1)>=0;
            % direction (depends on stable or unstable manifold)
            if Tmax < 0 % stable manifold
                dir1 = 1;
            else        % unstable manifold
                dir1 = -1;
            end
        end
        % if there is a second Poincare section provided
        if numarg == 2
            if strcmp(PS2,'U1')
                val2 = y(1,1)<0 && y(2,1)>=0;
                % direction (depends on stable or unstable manifold)
                if Tmax < 0 % stable manifold
                    dir2 = 1;
                else        % unstable manifold
                    dir2 = -1;
                end
            elseif strcmp(PS2,'U2')
                val2 = y(1,1)<=1-mu && y(2,1)<0;
                % direction (depends on stable or unstable manifold)
                if Tmax < 0 % stable manifold
                    dir2 = -1;
                else        % unstable manifold
                    dir2 = 1;
                end
            elseif strcmp(PS2,'U3')
                val2 = y(1,1)>=1-mu && y(2,1)>0; 
                % direction (depends on stable or unstable manifold)
                if Tmax < 0 % stable manifold
                    dir2 = 0;
                else        % unstable manifold
                    dir2 = -1;
                end
            elseif strcmp(PS2,'U4')
                val2 = y(1,1)-1<0 && y(2,1)>=0;
                % direction (depends on stable or unstable manifold)
                if Tmax < 0 % stable manifold
                    dir2 = -1;
                else        % unstable manifold
                    dir2 = 1;
                end
            end
        end
        % construct event condition assignments
        if numarg == 1
            value = val1;
            isterminal = 1;
            direction = dir1;
        elseif numarg == 2
            value = horzcat(val1,val2);
            isterminal = [1, 1];
            direction = horzcat(dir1,dir2);
        end     
    end
    
    % CR3BP Equation of Motion
	function dydt = ode45_CR3BP(t,X)
        % unpack state
        x = X(1,1);  y = X(2,1);  z = X(3,1);
        vx = X(4,1); vy = X(5,1); vz = X(6,1);
        % radii in the CR3BP
        r1 = sqrt( (x+mu)^2 + y^2 + z^2 );
        r2 = sqrt( (x-1+mu)^2 + y^2 + z^2 );
        % velocities
        dydt(1,1) = vx;
        dydt(2,1) = vy;
        dydt(3,1) = vz;
        % accelerations
        dydt(4,1) = 2*vy + x - ((1-mu)/r1^3)*(mu+x) + (mu/r2^3)*(1-mu-x);
        dydt(5,1) = -2*vx + y - ((1-mu)/r1^3)*y - (mu/r2^3)*y;
        dydt(6,1) = -((1-mu)/r1^3)*z - (mu/r2^3)*z;
    end

end
