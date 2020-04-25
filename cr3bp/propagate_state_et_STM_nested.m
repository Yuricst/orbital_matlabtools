function [rr,vv,time,stmcell] = ...
    propagate_state_et_STM_nested(mu_cr3bp,X0_input,Tf,nsteps,varargin)
% Nested function to propagate X0 in the CR3BP using ode45
% Yuri Shimane, 2020/04/07
% ======================================================= %
% FORMAT 
%   [rr,vv,time,stmcell] = ...
%      propagate_state_et_STM_nested(mu_cr3bp,X0,Tf,nsteps,...
%                                    'on',1e-5,1e-8)
% INPUT 
%   mu_cr3bp : CR3BP mass parameter
%   X0_input : 6x1 state-vector
%   Tf : final time (non-dimensional)
%   nsteps : number of steps 
%   (varargin includes)
%       event_option : 'off' or 'on'
%       relTol : relative error tolerance for ode45 option
%       absTol : absolute error tolerance for ode45 option
% OUTPUT
%   rr : N by 3 array of positions x,y,z (N is nsteps if 'off' Tf 
%           or resulting number of steps for 'on')
%   vv : N by 3 array of velocities xdot,ydot,zdot
%   time : N by 1 array of time (output from ode45)
%   stmcell : N by 1 cells each containing 6x6 STM
% ======================================================= %

% unpack varargin
if nargin == 4
    event_option = 'off';
    relTol = 1e-5;
    absTol = 1e-8;
elseif nargin == 5
    event_option = varargin{1};
    relTol = 1e-5;
    absTol = 1e-8;
elseif nargin > 5
    event_option = varargin{1};
    if nargin >= 6
        relTol = varargin{2};
        if nargin == 7
            absTol = varargin{3};
        end
    end
end

% convert X0 to double if it is a cell
if iscell(X0_input)==1
    X0init = zeros(6,1);
    for i = 1:6
        X0init(i,1) = X0_input{i,1};
    end
else
    X0init = X0_input;
end

% concatenate initial STM matrix to ode variable X
X0 = zeros(6,1); 
X0(1:6,1) = X0init;
stm_init = eye(6);
for row_init = 1:6
    for col_init = 1:6
        X0(6+col_init+(row_init-1)*6, 1) = stm_init(row_init,col_init);
    end
end

% declare global variable
global mu;
mu = mu_cr3bp;   % FIXME...?

% set-up solution space for stm
if strcmp(event_option, 'off')
    stmcell = cell(1, 1);
elseif strcmp(event_option, 'on')
    stmcell = cell(1, 1);
end
% initial STM
stm0 = eye(6);
stmcell{1,1} = stm0;
rr(1,:) = [X0(1,1), X0(2,1), X0(3,1)];
vv(1,:) = [X0(4,1), X0(5,1), X0(6,1)];
% create time array
time_arr = linspace(0, Tf, nsteps);
% run ode45
if strcmp(event_option, 'off')
    options = odeset('InitialStep',nsteps,...
        'RelTol',relTol,'AbsTol',absTol);
    [time, dynmat] = ode45(@rhs_cr3bp, time_arr, X0, options);
elseif strcmp(event_option, 'on')
    options = odeset('InitialStep',nsteps,'Events',@event_XZplane,...
        'RelTol',relTol,'AbsTol',absTol);
    [time, dynmat] = ode45(@rhs_cr3bp, [0 1.25*Tf], X0, options);
                                    % Tf is amplified to ensure condition
                                    % takes effect
end
% set-up solution space for rr and vv
rr = zeros(length(time), 3);
vv = zeros(length(time), 3);
% decompose dynmat to rr, vv, time, stmcell
rr(:,1) = dynmat(:,1);
rr(:,2) = dynmat(:,2);
rr(:,3) = dynmat(:,3);
vv(:,1) = dynmat(:,4);
vv(:,2) = dynmat(:,5);
vv(:,3) = dynmat(:,6);
for j = 1:length(time)
    stmmat = zeros(6,6);
    for myrow = 1:6
        for mycol = 1:6
            stmmat(myrow,mycol) = dynmat(j, 6+mycol+(myrow-1)*6);
        end
    end
    stmcell{j,1} = stmmat;
end

% ============== Nested Functions ============== %
    % Event function for stopping criteria with 'on' Tf option
    function [value,isterminal,direction] = event_XZplane(t, X)
        % when the y-value changes sign (=0)
        value = X(2,1);  % when value = 0
        isterminal = 1;
        direction = 0;
    end
    
    % CR3BP Equation of Motion
    function Xdot = rhs_cr3bp(t, X)
        % Function evaluates derivative of state
        % decompose state and STM
        x = X(1,1);  y = X(2,1);  z = X(3,1);
        vx = X(4,1); vy = X(5,1); vz = X(6,1);
        stm = zeros(6,6);
        for row = 1:6
            for col = 1:6
                stm(row,col) = X(6+col+(row-1)*6, 1);
            end
        end
        % calculate radii
        r1 = sqrt( (x+mu)^2 + y^2 + z^2 );
        r2 = sqrt( (x-1+mu)^2 + y^2 + z^2 );
        % --- STATE DERIVATIVE --- %
        % initialize double
        Xdot = zeros(42,1);
        % position-state derivative
        Xdot(1,1) = vx;
        Xdot(2,1) = vy;
        Xdot(3,1) = vz;
        % velocity-state derivative
        Xdot(4,1) = 2*vy + x - ((1-mu)/r1^3)*(mu+x) + (mu/r2^3)*(1-mu-x);
        Xdot(5,1) = -2*vx + y - ((1-mu)/r1^3)*y - (mu/r2^3)*y;
        Xdot(6,1) = -((1-mu)/r1^3)*z - (mu/r2^3)*z;
        % --- A-MATRIX --- %
        % initialize A-matrix
        A11 = zeros(3,3);
        A12 = eye(3);
        A22 = 2*[0  1 0;
                 -1 0 0;
                 0  0 0];
        % construct U matrix
        Uxx = 1 - (1-mu)/r1^3 - mu/r2^3 + 3*((x+mu)^2 *(1-mu)/r1^5 ...
                    + (x+mu-1)^2*mu/r2^5);
        Uyy = 1 - (1-mu)/r1^3 - mu/r2^3 + 3*y^2*((1-mu)/r1^5 + mu/r2^5);
        Uzz = -(1-mu)/r1^3 - mu/r2^3 + 3*z^2*((1-mu)/r1^5 + mu/r2^5);
        Uxy = 3*y*((x+mu)*(1-mu)/r1^5 + (x+mu-1)*mu/r2^5);
        Uxz = 3*z*((x+mu)*(1-mu)/r1^5 + (x+mu-1)*mu/r2^5);
        Uyz = 3*y*z*((1-mu)/r1^5 + mu/r2^5);
        Uderiv = [Uxx, Uxy, Uxz;
                  Uxy, Uyy, Uyz;
                  Uxz, Uyz, Uzz];
        % update A-matrix
        A = [A11,    A12;
             Uderiv, A22];
        % differential relation
        stmdot = A * stm;
        % store elements of A in Xdot(7,1) ~ onwards
        for row = 1:6
            for col = 1:6
                Xdot(6+col+(row-1)*6, 1) = stmdot(row,col);
            end
        end
    end
% ============================================== %
end