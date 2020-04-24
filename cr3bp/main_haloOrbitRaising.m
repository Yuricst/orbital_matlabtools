%% main - optimal control same-libration-point halo orbit raising
% main function to setup and solve optimal control problem for
% orbit raising a halo orbit about a libration point 
% Yuri Shimane, 2020/04/17

% house keep
clear; close all; clc;
tic;

% define CR3BP parameter
defineParam_EarthMoon;

%% Halo construction (orbit 1)
% ... choose in-plane phase angle
phi = pi;
% ... choose collinear Lagrangan point (1,2, or 3)
lp = 1;
% ... choose northern or southern family
northsouth = 3;  % n = 1 (southern) or 3 (northern)
% ... choose out-of-plane altitude of halo-orbit [km]
Az_km = [100, 110]; %linspace(0,2000,10)';
% analytical 3rd order solution
fprintf('Create analytical 3rd order solution...\n');
% construct Halo orbit
[X0_halo1,Thalo1,~,~] = ...
    constructhaloX0(mu,Lstar,phi,lp,northsouth,Az_km(1));
[X0_halo2,Thalo2,~,~] = ...
    constructhaloX0(mu,Lstar,phi,lp,northsouth,Az_km(2));
% propagate final result for plotting
[rr_halo1,vv_halo1,time_halo1,stmcell_halo1] = ...
    propagate_state_et_STM_nested(mu,X0_halo1,1*Thalo1,4000);
[rr_halo2,vv_halo2,time_halo2,stmcell_halo2] = ...
    propagate_state_et_STM_nested(mu,X0_halo2,1*Thalo2,4000);
% plot result
figure(98)
plot3(LP(lp,1),LP(lp,2),LP(lp,3),'ok');
hold on;
plot3(rr_halo1(:,1), rr_halo1(:,2), rr_halo1(:,3),'b');
hold on
plot3(rr_halo2(:,1), rr_halo2(:,2), rr_halo2(:,3),'r');
grid on;
xlabel('x'); ylabel('y'); zlabel('z');

%% sandbox... play around with angular momentum
% compute angular momentum vector
hvect1 = zeros(3,length(time_halo1));
hnorm1 = zeros(length(time_halo1),1);
for i = 1:length(time_halo1)
    rtmp = rr_halo1(i,:)';
    vtmp = vv_halo1(i,:)';
    htmp = cross(rtmp, vtmp);
    hvect1(:,i) = htmp;  
    hnorm1(i,1) = norm(hvect1(:,i));
end

figure(2)
plot(time_halo1,hnorm1);
grid on;


%% Setup fmincon
opts = optimoptions('fmincon','Display','iter',...
        'MaxFunctionEvaluations',3e5,...
        'MaxIterations',4e3,'ConstraintTolerance',1.0000e-05,...
        'FiniteDifferenceType','forward','FiniteDifferenceStepSize',1e-8);
    
% number of steps
nsteps = 100;
% initial condition
u0 = [ 0.8 * ones(nsteps,1);   % initial guess for u1 -- x-direction
       0.1 * ones(nsteps,1);   % initial guess for u2 -- y-direction
       0 * ones(nsteps,1);   % initial guess for u3 -- z-direction
       3*Thalo2];            % initial guess for time of flight
% initial halo orbit state
X0 = X0_halo1;
% target halo orbit state
Xf = X0_halo2;
% max time of flight
Tfmax = 4*Thalo2;
% initial mass of s/c
mass0 = 500;  % [kg]
% Thruster: RIT-10 EVO
% Isp
Isp = 3400;     % [sec]
% max F
Fmax = 15e-3; % [N]
toc;
disp('Running fmincon...');
% run nested optimization function
[u,fval,exitflag,output] = ...
    runoptim_raiseOrb(mu,u0,mass0,X0,Xf,Isp,Fmax,Tfmax,opts,nsteps,Lstar,Tstar);

toc;

%% Plot result
% propagate result with ode45() for plotting
[time,dynmat] = propagate_cr3bp_thrust(mu,u,mass0,X0,Isp,Fmax,4000);
figure(3)
plot3(dynmat(:,1),dynmat(:,2),dynmat(:,3));
grid on;

figure(4)
plot(time,dynmat(:,7));
grid on;
xlabel('Time'); ylabel('mass [kg]');
% ... check mass dynamics?















