%% testing bed

% house keeping
clear; close all; clc;

nsteps = 100;
u = ones(200,1);
time = linspace(0,3.32,nsteps);

%% ============== %

nsteps = length(time);
dt = time(1,2) - time(1,1);
% handle for rhs expressions of dynamical system
rhs = @rhs_orbit_raising;

% initialize cells to store
rcell      = cell(nsteps,1);
thetacell  = cell(nsteps,1);
vrcell     = cell(nsteps,1);
vthetacell = cell(nsteps,1);
% initial conditions in cell
rcell{1,1}      = 1;
thetacell{1,1}  = 0;
vrcell{1,1}     = 0;
vthetacell{1,1} = 1;
% initialize array to store real value
r      = zeros(nsteps,1);
theta  = zeros(nsteps,1);
vr     = zeros(nsteps,1);
vtheta = zeros(nsteps,1);

% MCX integration breakdown
mcx_discretization = 2*nsteps/5;  % -- divide total steps by 5 (5 mcx dir at a time)
% MCX step size
h = 1e-10;
% initialize and populate cell holding real control values
ucell = cell(2*nsteps,1);
for i = 1:2*nsteps
    ucell{i,1} = u(i,1);
end

tic;

% integration is repeated to obtain all necessary gradients
for j = 1:mcx_discretization  % -- runs through from 1:40 (1:5:200)
    % prepare MCX arrays: input to function multicomplex()
    u_p = zeros(5,2^5);
    for k = 1:5
        % populate first element with real value of u
        u_p(k,1) = u((j-1)*5 + k, 1);
        % add perturbation in k-th direction
        u_p(k,1+2^(k-1)) = h;
    end
    
    % replace real u with MCX u
    for i = 1:2*nsteps
        if i == (j-1)*5 + 1
            ucell{i,1} = multicomplex(u_p(1,:));
        elseif i == (j-1)*5 + 2
            ucell{i,1} = multicomplex(u_p(2,:));
        elseif i == (j-1)*5 + 3
            ucell{i,1} = multicomplex(u_p(3,:));
        elseif i == (j-1)*5 + 4
            ucell{i,1} = multicomplex(u_p(4,:));
        elseif i == (j-1)*5 + 5
            ucell{i,1} = multicomplex(u_p(5,:));
        end
    end
    
    % integrate dynamics (with mcx vector when relevant)
    for i = 1:nsteps-1
        % compute rhs
        [rdot,thetadot,vrdot,vthetadot] = rhs(rcell{i,1},thetacell{i,1},vrcell{i,1},...
            vthetacell{i,1},ucell{i,1},ucell{nsteps+i,1},time(i));
        % update with time-step (EULER method)
        rcell{i+1,1}      = rcell{i,1}      + dt*rdot;
        thetacell{i+1,1}  = thetacell{i,1}  + dt*thetadot;
        vrcell{i+1,1}     = vrcell{i,1}     + dt*vrdot;
        vthetacell{i+1,1} = vthetacell{i,1} + dt*vthetadot;
        
    end
    
    % extract sensitivities (gradients) and store
    tmp = 0;
    
    % replace MCX u back with real u before next integration
    for i = 1:2*nsteps
        if i == (j-1)*5 + 1
            ucell{i,1} = u(i,1);
        elseif i == (j-1)*5 + 2
            ucell{i,1} = u(i,1);
        elseif i == (j-1)*5 + 3
            ucell{i,1} = u(i,1);
        elseif i == (j-1)*5 + 4
            ucell{i,1} = u(i,1);
        elseif i == (j-1)*5 + 5
            ucell{i,1} = u(i,1);
        end
    end
end

toc;

% dynmat = [r, theta, vr, vtheta];

gradf = 0;


