%% main orbit raising - nested
% Yuri Shimane, 2020/02/20
% run orbit raising problem with nested functions fed to fmincon
% ====================================== %

% house keeping
clear; close all; clc;
tic;

% optimizer options - decide whether to use MCX
mymcxop = false; % true OR false 

if mymcxop == false
    % 1. default, fmincon uses its own gradient calculation  
    opts = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',1.5e5,... % <1.5e5 gives optimal solution
        'MaxIterations',4e3,'ConstraintTolerance',1.0000e-03,...
        'FiniteDifferenceType','forward','FiniteDifferenceStepSize',1e-8);
elseif mymcxop == true
    % 2. feeding mcx-computed gradients & sensitivities
    disp('Using MCX for gradient evaluation;')
    % >>> follows syntax given by: 
    % https://uk.mathworks.com/help/optim/ug/fmincon-interior-point-algorithm-with-analytic-hessian.html
    opts = optimoptions('fmincon','Display','iter-detailed','MaxFunctionEvaluations',100,...
        'MaxIterations',2,'ConstraintTolerance',1.0000e-06,...  % -- allow only one iter for now!
        'CheckGradients',false,...
        'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',false);
end

% number of steps
nsteps = 50;
% initial condition
u0 = [ 0 * ones(nsteps,1);...   % initial guess for u1 -- radial
       1 * ones(nsteps,1)];     % initial guess for u2 -- theta
% run nested optimization function
[u,fval,exitflag,output,c,ceq] = runoptim(u0,opts,nsteps,mymcxop);
toc;
% Iter displays: https://uk.mathworks.com/help/optim/ug/iterative-display.html#f92519
% Feasibility > Maximum constraint violation, where satisfied inequality constraints count as 0
% First-order optim > First-order optimality measure (should be 0)

% check for constraints
disp('Done!')
fprintf('Final altitude: %6.4f\n',-fval)
fprintf('vr req @tf: %6.4f, vtheta req @tf: %6.4f\n',ceq(1,5),ceq(1,6))


%% Plot result
% deconstruct result
u1 = u(1:nsteps, 1);
u2 = u(nsteps+1:length(u), 1);
% time
tf = 3.32;
% discretizedtime array
time = linspace(0,tf,nsteps);
% recreate orbit
dynmat_sol = dynamics_orbit_raising(time,u);
r = dynmat_sol(:,1);
theta = dynmat_sol(:,2);
vr = dynmat_sol(:,3);
vtheta = dynmat_sol(:,4);
% create control vector in cartesian
ux = u1.*cos(theta) - u2.*sin(theta);
uy = u1.*sin(theta) + u2.*cos(theta);

% plot results
figure(10)
hold on
plot(time,r,'-xb');
plot(time,theta,'-xr');
plot(time,vr,'-xc');
plot(time,vtheta,'-xm');
legend('r','\theta','v_r','v_{\theta}','location','northwest')
xlabel('time');
grid on; grid minor;

figure(11)
u_total = u1.^2 + u2.^2;  % -- total control
plot(time,u1,'-x');
hold on
plot(time,u2,'-x');
hold on
plot(time,u_total,'-x');
ylim([-1 1])
legend('u1: ur','u2: utheta','u1^2+u2^2')
xlabel('time');
grid on; grid minor;

figure(12)
subplot(2,1,1)
hold on
plot(time,r,'-x')
plot(time,theta,'-x')
legend('r','theta')
xlabel('time');
grid on; grid minor;

subplot(2,1,2)
hold on
plot(time,vr,'-x')
plot(time,vtheta,'-x')
legend('vr','vtheta')
xlabel('time');
grid on; grid minor;

% orbit plot
for i = 1:nsteps
    xy_traj(i,:) = r(i,1) * [cos(theta(i,1)), sin(theta(i,1))];
end
% initial and final positions
theta_ref = linspace(0,2*pi,100)';
for i = 1:100
    xy_initial(i,:) = r(1,1) * [cos(theta_ref(i,1)), sin(theta_ref(i,1))];
    xy_final(i,:) = r(length(r),1) * [cos(theta_ref(i,1)), sin(theta_ref(i,1))];
end

figure(13)
hold on
% plot transfer trajectory
plot(xy_traj(:,1), xy_traj(:,2),'-x');
plot(xy_traj(1,1), xy_traj(1,2),'*');  % marker at initial position
plot(xy_traj(nsteps,1), xy_traj(nsteps,2),'^'); % marker at final position
% initial and final trajectory
plot(xy_initial(:,1), xy_initial(:,2),'-.k');
plot(xy_final(:,1), xy_final(:,2),'--k');
% control vector
quiver(xy_traj(:,1), xy_traj(:,2),ux, uy)
caption_finalorbit = strcat('final circular @r = ', num2str(-fval));
legend('transfer','start','end','initial circular',caption_finalorbit,'accel. vector',...
    'location','southwest')
axis equal
grid on; grid minor;
xlabel('x'); ylabel('y');

%% analyze evolution of (some) orbital elements?
% >>> particularly eccentricity (=semi-major axis?)

%% dimensionalize solution
% physical parameters
AU = 149598870.700; % 1 AU in km
mu = 132712440018; % sun gravitational parameter
% non-dimensionalization units
r0_phys = 1*AU;
v0_phys = sqrt(mu/r0_phys);
% positions
xy_traj_phys = r0_phys * xy_traj;
xy_initial_phys = r0_phys * xy_initial;
xy_final_phys = r0_phys * xy_final;
r_phys = r0_phys * r;
% time
time_phys = time * sqrt(r0_phys/mu);
% velocities
vr_phys = v0_phys * vr;
vtheta_phys = v0_phys * vtheta;

% compute initial and final orbits velocities
vcirc_0  = sqrt(mu/r_phys(1,1));
vcirc_fn = sqrt(mu/r_phys(nsteps,1));

% figure(21)
% hold on
% plot(time/(60*60*24), vr_phys);
% plot(time/(60*60*24), vtheta_phys);
% yline(vcirc_0,'-.k');
% yline(vcirc_fn,'--k');
% legend('transfer v_r','transfer v_{\theta}','initial orbit v_{circ}','final orbit v_{circ}',...
%     'location','west');
% xlabel('time [days]'); ylabel('v [km/s]');
% grid on; grid minor;
% 
% figure(22)
% hold on
% plot(xy_traj_phys(:,1),xy_traj_phys(:,2),'-x')
% plot(xy_initial_phys(:,1),xy_initial_phys(:,2),'-.k')
% plot(xy_final_phys(:,1),xy_final_phys(:,2),'--k')
% grid on; grid minor;
% axis equal
% xlabel('x [km]'); ylabel('y [km]');


