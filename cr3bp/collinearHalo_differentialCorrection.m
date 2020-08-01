function [X0_halo,Tfull,varargout] = ...
    collinearHalo_differentialCorrection(mu,X0_input,Thalf,residualtol,varargin)
% Function applies differential correction to state X0 to generate 
% periodic halo orbit about collinear libratin point
% Yuri Shimane, 2020/04/08
% ======================================================= %
% INPUT
%   mu : cr3bp mass parameter
%   X0_DC : 6x1 array of guessed initial state of halo
%   Thalf : initial guess for half-period of halo (Thalf ~ T/2)
%   residualtol : tolerance on residuals on dvx and dvz at X at Thalf
%   nsteps : (default: 4000) approximate number of steps to integrate with ode45
%   reltol : (default: 1e-5) relative tolerance for ode45
%   abstol : (default: 1e-7) absolute tolerance for ode45
%   DC_lim : (default: 6) # of DC iteration allowed
% OUTPUT
%   X0_halo : 6x1 array of initial state for Halo
%   Tfull : period
%   rr : (optional) N x 3 position array of corrected X0 propagated until Thalf
%   vv : (optional) N x 3 velocity array of corrected X0 propagated until Thalf
%   time : (optional) N x 1 array of time
%   stmcell : (optional) N x 1 cell containing STM from t=0 to each time step
%   residual_log : (optional) 2 x (# of DC iteration) log of residuals dvx and dvz
% ======================================================= %

% unpack varargin
if nargin == 4
    nsteps = 4000;
    reltol = 1e-5;
    abstol = 1e-7;
    DC_lim = 12;
elseif nargin == 5
    nsteps = varargin{1};
    reltol = 1e-5;
    abstol = 1e-7;
    DC_lim = 12;
elseif nargin == 6
    nsteps = varargin{1};
    reltol = varargin{2};
    abstol = 1e-7;
    DC_lim = 12;
elseif nargin == 7
    nsteps = varargin{1};
    reltol = varargin{2};
    abstol = varargin{3};
    DC_lim = 12;
elseif nargin == 8
    nsteps = varargin{1};
    reltol = varargin{2};
    abstol = varargin{3};
    DC_lim = varargin{4};
end

% convert X0 to double if it is a cell
if iscell(X0_input)==1
    X0_DC = zeros(6,1);
    for i = 1:6
        X0_DC(i,1) = X0_input{i,1};
    end
else
    X0_DC = X0_input;
end

% initialize DC residue
residue_DC = 10 * residualtol;
residual_log = zeros(2,1);
% setup counter
count = 1;
while residue_DC > residualtol
    % print message
    fprintf('Differential correction iteration %d\n',count);
    % propagate
    [rr,vv,time,stmcell] = ...
        propagate_state_et_STM_nested(mu,X0_DC,2*Thalf,nsteps,'on',reltol,abstol);
    % check final residuals dvx and dvz
    residual_vec = [vv(length(time), 1);
                    vv(length(time), 3)];
    % log residual
    residual_log(:,count) = residual_vec(:,1);
    % STM after half-period
    stm_T2 = stmcell{length(time),1};
    % final time-step
    dt_tf = time(length(time),1) - time(length(time)-1,1);
    % compute xdotdot, zdotdot, ydot at half-period (backward euler)
    ydot_tf = vv(length(time),2);
    xdotdot = ( vv(length(time),1) - vv(length(time)-1,1) )/dt_tf;
    zdotdot = ( vv(length(time),3) - vv(length(time)-1,3) )/dt_tf;
    % Construct differential correction map [x0; ydot0] -> [xdot; zdot]
    D1 = [stm_T2(4,1) stm_T2(4,5);
          stm_T2(6,1) stm_T2(6,5)];
    D2 = -(1/ydot_tf) * [xdotdot; zdotdot] * [stm_T2(2,1), stm_T2(2,5)];
    B_map0f = D1 + D2;
    % compute residue at T = 0
    residue0 = B_map0f \ residual_vec;
    % print result
    fprintf('Residue dvx:  %2.8s, dvz: %2.8s \n',residual_vec(1,1),residual_vec(2,1))
    fprintf('Correction x: %2.8s, vy:  %2.8s \n',residue0(1,1),residue0(2,1));
    % modify initial condition
    X0_DC(1,1) = X0_DC(1,1) - residue0(1,1);  % update x
    X0_DC(5,1) = X0_DC(5,1) - residue0(2,1);  % update vy
    % update half-period length for next iteration
    Thalf = time(end);
    % store larger of the two residual
    if abs(residual_vec(1,1)) > abs(residual_vec(2,1))
        residue_DC = abs(residual_vec(1,1));
    else
        residue_DC = abs(residual_vec(2,1));
    end
    % display error if number of DC iteration exceeds limit
    if count > DC_lim
        error('Error: exceeded allowable DC iteration set at %d\n',count);
    end
    % update counter
    count = count + 1;
end

% extract final iteration initial state
X0_halo = [rr(1,1);
           rr(1,2);
           rr(1,3);
           vv(1,1);
           vv(1,2);
           vv(1,3)];
       
% prepare variable output arguments
Tfull = 2*time(end);  % full period of Halo
varargout{1} = rr;
varargout{2} = vv;
varargout{3} = time;
varargout{4} = stm_T2;
varargout{5} = stmcell;
varargout{6} = residual_log;

end
