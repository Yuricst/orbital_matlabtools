%% main script for EP-powered halo
% Yuri Shimane, 2020/05/02
% Script plots EP-powered manifold for prelim. investigation
% ------------------------------------------- %

% house keep
clear; close all; clc;

% define CR3BP parameter
defineParam_EarthMoon;

%% Construct Halo orbit
% choose Halo parameters
phi = 0;  % in-plane phase angle
lp = 1;  % choose collinear Lagrangian point (1, 2, or 3)
northsouth = 1;  % n = 1 (southern) or 3 (northern)
Az_km = 26000;  % amplitude [km] 0.01 / 100 / 200 / 8000 / 26000
% construct Halo orbit
[X0_halo,Thalo,stm_T2_halo,~] = ...
    constructhaloX0(mu,Lstar,phi,lp,northsouth,Az_km);
% propagate final result for plotting
[rr_halo,vv_halo,time_halo,stmcell_halo1] = ...
    propagate_state_et_STM_nested(mu,X0_halo,1*Thalo,8000);
% Jacobi constant of Halo orbit
[halo_C, ~] = jacobiConst(mu, rr_halo, vv_halo);

% plot result
figure(1)
plot3(LP(lp,1),LP(lp,2),LP(lp,3),'xk');
hold on;
plot3(rr_halo(:,1), rr_halo(:,2), rr_halo(:,3),'b');
grid on;
xlabel('x'); ylabel('y'); zlabel('z');
title(['Earth-Moon Halo, Az = ',num2str(Az_km),' [km]']);

% plot Halo with primaries and L point
figure(101)
hold on
plot(-mu,0,'ob');
plot(1-mu,0,'ok');
plot(LP(lp,1),LP(lp,2),'xk');
plot(rr_halo(:,1), rr_halo(:,2),'-b','LineWidth',1.5);

[Xsphere, Ysphere, Zsphere] = sphere(20);

figure(102)
mesh(6378*Xsphere/Lstar-mu, 6378*Ysphere/Lstar, 6378*Zsphere/Lstar, 'edgecolor', 'b');
hold on
mesh(1737*Xsphere/Lstar+1-mu, 1737*Ysphere/Lstar, 1737*Zsphere/Lstar, 'edgecolor', 'k');
plot3(LP(lp,1),LP(lp,2),LP(lp,3),'xk');
plot3(rr_halo(:,1), rr_halo(:,2), rr_halo(:,3),'-b','LineWidth',2);


%% Define thrust and spacecraft mass
% spacecraft
m_sc0 = 500; % spacecraft mass [kg]
mass0_nd = m_sc0 / m_sc0;  % spacecraft mass [non-dim]
% thruster
P = 435; % [W]
Isp = 3400;  % [sec]
eta = 0.52;  % efficiency (RIT-10 has 0.52)
Thrust = 2*eta*P/(9.81*Isp);
fprintf('RIT-10 EVO thrust: %4.6f [N]\n',Thrust);
% final mass
mdot = -2*eta*P/(9.81*Isp)^2;
mdot_nd = mdot / (mstar/Tstar);
f = Thrust / (mass0_nd*Lstar/Tstar^2); % should be ~e-2

f = 6.95e-2;  % FORCE-FIX to deep-space one
mdot_nd = -f*Lstar/(Isp*9.81*Tstar);

% >>> https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Conferences/2018_IAC_PriCoxHowFolGre.pdf
% >>> https://engineering.purdue.edu/people/andrew.d.cox.2/pubs/2019-08_AAS_CoxHowFol.pdf


%% Construct manifolds
disp('Constructing manifolds...')
% construct monodromy matrix
[M1,eigval1,eigvec1] = constructMonodromy(stm_T2_halo);
% stable/unstable vector (normalized)
Yu = eigvec1(:,1)/norm(eigvec1(:,1));
Ys = eigvec1(:,2)/norm(eigvec1(:,2));
% perturbation of state
epsobj = 1e-4;  % eps objective
% eps factor for linear approximated initial guess
eigvec_average = (abs(Yu(1,1)) + abs(Yu(3,1)))/2;
eps = epsobj/eigvec_average;
% decide number of manifolds to make
num_manif = 30;
% create manifolds
for i = 1:num_manif
    % obtain num_manif points along the halo
    if i == 1
        manif(i).X0index = 1;  % ugly fix...
    else
        manif(i).X0index = (i-1)*floor(length(time_halo)/num_manif);
    end
    manif(i).X0 = [rr_halo(manif(i).X0index, 1);
                   rr_halo(manif(i).X0index, 2);
                   rr_halo(manif(i).X0index, 3);
                   vv_halo(manif(i).X0index, 1);
                   vv_halo(manif(i).X0index, 2);
                   vv_halo(manif(i).X0index, 3)];
    % transcribe eigenvectors using STM
    manif(i).Yu = stmcell_halo1{manif(i).X0index,1}*Yu/norm(stmcell_halo1{manif(i).X0index,1}*Yu);
    manif(i).Ys = stmcell_halo1{manif(i).X0index,1}*Ys/norm(stmcell_halo1{manif(i).X0index,1}*Ys);
    % unstable linear perturbed states
    manif(i).X0u_p = manif(i).X0 + eps*manif(i).Yu;
    manif(i).X0u_m = manif(i).X0 - eps*manif(i).Yu;
    % stable linear perturbed states
    manif(i).X0s_p = manif(i).X0 + eps*manif(i).Ys;
    manif(i).X0s_m = manif(i).X0 - eps*manif(i).Ys;


    % ============ Regular Manifold ============ %
    % propagate stable manifolds (backward in time)
    [manif(i).rr_s_p, manif(i).vv_s_p, manif(i).time_s_p] = ...
        propagate_manifold_poincare(mu,manif(i).X0s_p,-3*Thalo,'U1','U3');
%     [manif(i).rr_s_m, manif(i).vv_s_m, manif(i).time_s_m] = ...
%         propagate_manifold_poincare(mu,manif(i).X0s_m,-4*Thalo1,'U1','U3');
    
    % ============ EP powered ============ %
    % time to turn Thrust on after propagation begins (abs value)
    t_thrustON = 0.4*Thalo;

    % EP powered manifold
    [manif(i).EPrr_s_p, manif(i).EPvv_s_p, ...
     manif(i).EPtime_s_p, manif(i).EPmass_s_p] = ...
        propagate_EPmanifold_poincare(mu,manif(i).X0s_p,-3*Thalo,...
                                        f,t_thrustON,mass0_nd,mdot_nd,'U1','U3');
%     [manif(i).EPrr_s_m, manif(i).EPvv_s_m, ...
%      manif(i).EPtime_s_m, manif(i).EPmass_s_m] = ...
%         propagate_EPmanifold_poincare(mu,manif(i).X0s_m,-0.1*Thalo1,...
%                                         Thrust_nd,mass0_nd,mdot_nd,'U1','U3');

    % ============ PLOTS ============ %
    % append to XY plot of manifolds
    figure(101)
    hold on
    % plot manifolds
    plot(manif(i).rr_s_p(:,1),manif(i).rr_s_p(:,2),'-g');
    plot(manif(i).EPrr_s_p(:,1),manif(i).EPrr_s_p(:,2),'-m');
    grid on;
    xlabel('x'); ylabel('y');
    
    % 3D plot of manifolds
    figure(102)
    hold on
    % plot manifolds
    plot3(manif(i).rr_s_p(:,1),manif(i).rr_s_p(:,2),manif(i).rr_s_p(:,3),'-g');
    plot3(manif(i).EPrr_s_p(:,1),manif(i).EPrr_s_p(:,2),manif(i).EPrr_s_p(:,3),'-m');
    grid on;
    xlabel('x'); ylabel('y');

    
end

figure(101)
legend('Earth','Moon','L1','Halo','nominal manif.','powered manif.')
figure(102)
legend('Earth','Moon','L1','Halo','nominal manif.','powered manif.')


%% Plot
figure(31)
subplot(2,1,1)
hold on;
% nominal manifold
plot(manif(1).time_s_p, manif(1).rr_s_p(:,1), '-b');
plot(manif(1).time_s_p, manif(1).rr_s_p(:,2), '-r');
plot(manif(1).time_s_p, manif(1).rr_s_p(:,3), '-m');
% EP manifold
plot(manif(1).EPtime_s_p, manif(1).EPrr_s_p(:,1), '--b');
plot(manif(1).EPtime_s_p, manif(1).EPrr_s_p(:,2), '--r');
plot(manif(1).EPtime_s_p, manif(1).EPrr_s_p(:,3), '--m');
% Thruster ON/OFF time
xline(-t_thrustON,'--k');
grid on;

subplot(2,1,2)
hold on;
% nominal manifold
plot(manif(1).time_s_p, manif(1).vv_s_p(:,1), '-b');
plot(manif(1).time_s_p, manif(1).vv_s_p(:,2), '-r');
plot(manif(1).time_s_p, manif(1).vv_s_p(:,3), '-m');
% EP manifold
plot(manif(1).EPtime_s_p, manif(1).EPvv_s_p(:,1), '--b');
plot(manif(1).EPtime_s_p, manif(1).EPvv_s_p(:,2), '--r');
plot(manif(1).EPtime_s_p, manif(1).EPvv_s_p(:,3), '--m');
% Thruster ON/OFF time
xline(-t_thrustON,'--k');
grid on;













