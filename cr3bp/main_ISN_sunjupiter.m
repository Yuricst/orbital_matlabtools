%% Sun-Jupiter for ITN
% Yuri Shimane, 2020/05/23

% house keeping
clear; close all; clc;

% define CR3BP parameter
defineParam_SunJupiter;

%% Halo construction
% ... choose in-plane phase angle
phi = 0;
% ... choose collinear Lagrangan point (1,2, or 3)
lp = 1;
% ... choose northern or southern family
northsouth = 1;  % n = 1 (southern) or 3 (northern)
% ... choose out-of-plane altitude of halo-orbit [km]
Az_km = 52000; %43800 is closest approach [Rausch 2005]
% analytical 3rd order solution
fprintf('Create analytical 3rd order solution...\n');
% construct Halo orbit
[X0_halo1,Thalo1,stm_T2_halo1,~] = ...
    constructhaloX0(mu,Lstar,phi,lp,northsouth,Az_km);
% propagate final result for plotting
[rr_halo1,vv_halo1,time_halo1,stmcell_halo1] = ...
    propagate_state_et_STM_nested(mu,X0_halo1,1*Thalo1,8000);
% Jacobi constant of Halo orbit
[halo_C, ~] = jacobiConst(mu, rr_halo1, vv_halo1);

% plot result
figure(1)
plot3(LP(lp,1),LP(lp,2),LP(lp,3),'ok');
hold on;
plot3(rr_halo1(:,1), rr_halo1(:,2), rr_halo1(:,3),'b');
grid on;
xlabel('x'); ylabel('y'); zlabel('z');
title(['Earth-Moon Halo, Az = ',num2str(Az_km),' [km]']);


%% Construct manifolds
disp('Constructing manifolds...')
% construct monodromy matrix
[M1,eigval1,eigvec1] = constructMonodromy(stm_T2_halo1);
% stable/unstable vector (normalized)
Yu = eigvec1(:,1)/norm(eigvec1(:,1));
Ys = eigvec1(:,2)/norm(eigvec1(:,2));   % NEED CHECKING!
% perturbation of state
epsobj = 1e-6;  % eps objective
% eps factor for linear approximated initial guess
eigvec_average = (abs(Yu(1,1)) + abs(Yu(3,1)))/2;
eps = epsobj/eigvec_average;
% decide number of manifolds to make
num_manif = 50;
% create manifolds
for i = 1:num_manif
    % obtain num_manif points along the halo
    if i == 1
        manif(i).X0index = 1;  % ugly fix...
    else
        manif(i).X0index = (i-1)*floor(length(time_halo1)/num_manif);
    end
    manif(i).X0 = [rr_halo1(manif(i).X0index, 1);
                   rr_halo1(manif(i).X0index, 2);
                   rr_halo1(manif(i).X0index, 3);
                   vv_halo1(manif(i).X0index, 1);
                   vv_halo1(manif(i).X0index, 2);
                   vv_halo1(manif(i).X0index, 3)];
    % transcribe eigenvectors using STM
    manif(i).Yu = stmcell_halo1{manif(i).X0index,1}*Yu/norm(stmcell_halo1{manif(i).X0index,1}*Yu);
    manif(i).Ys = stmcell_halo1{manif(i).X0index,1}*Ys/norm(stmcell_halo1{manif(i).X0index,1}*Ys);
    % unstable linear perturbed states
    manif(i).X0u_p = manif(i).X0 + eps*manif(i).Yu;
    manif(i).X0u_m = manif(i).X0 - eps*manif(i).Yu;
    % stable linear perturbed states
    manif(i).X0s_p = manif(i).X0 + eps*manif(i).Ys;
    manif(i).X0s_m = manif(i).X0 - eps*manif(i).Ys;
    % propagate stable manifolds (backward in time)
    [rr_maniftmp_s_p, vv_maniftmp_s_p, time_maniftmp_s_p] = ...
        propagate_manifold_poincare(mu,manif(i).X0s_p,-4*Thalo1,'U1','U3');
    [rr_maniftmp_s_m, vv_maniftmp_s_m, time_maniftmp_s_m] = ...
        propagate_manifold_poincare(mu,manif(i).X0s_m,-4*Thalo1,'U1','U3');
    
    % store
    manif(i).rr_s_p = rr_maniftmp_s_p;
    manif(i).vv_s_p = vv_maniftmp_s_p;
    manif(i).time_s_p = time_maniftmp_s_p;
    manif(i).rr_s_m = rr_maniftmp_s_m;
    manif(i).vv_s_m = vv_maniftmp_s_m;
    manif(i).time_s_m = time_maniftmp_s_m;
    % find closest approach of manifold to first primary
    for j = 1:length(rr_maniftmp_s_p)
        manif(i).rr_manif_norm_s_p(j) = norm(rr_maniftmp_s_p(j,:));
    end
%     for j = 1:length(rr_maniftmp_s_m)
%         manif(i).rr_manif_norm_s_m(j) = norm(rr_maniftmp_s_m(j,:));
%     end
    % (pls family)
    [manif(i).rr_min_s_p, manif(i).rr_minIndx_s_p] = min(manif(i).rr_manif_norm_s_p);
    % (min family)
%     [manif(i).rr_min_s_m, manif(i).rr_minIndx_s_m] = min(manif(i).rr_manif_norm_s_m);
    
    % angular momentum along the manifold
    % (pls family)
    [manif(i).h_s_p , manif(i).h_norm_s_p] = ...
        angularMomentum(rr_maniftmp_s_p, vv_maniftmp_s_p);
    % (min family)
    [manif(i).h_s_m , manif(i).h_norm_s_m] = ...
        angularMomentum(rr_maniftmp_s_m, vv_maniftmp_s_m);
    
    % inclination along the manifold
    manif(i).inc_arr = inclinationArray(manif(i).h_s_p);
    
    % Jacobi constant
    [manif(i).C_p, ~] = jacobiConst(mu, rr_maniftmp_s_p, vv_maniftmp_s_p);
    [manif(i).C_m, ~] = jacobiConst(mu, rr_maniftmp_s_m, vv_maniftmp_s_m);
    
    % plot of h_norm and C of each manifold
    figure(5)
    subplot(5,1,1)
    plot(manif(i).time_s_p, manif(i).h_norm_s_p, '-m');
    hold on;
    grid on;
    xlabel('time'); ylabel('h (manifold)');
    subplot(5,1,2)
    plot(manif(i).time_s_p, manif(i).C_p, '-m');
    hold on;
    grid on;
    xlabel('time'); ylabel('Jacobi Constant C');
    subplot(5,1,3)
    plot(manif(i).time_s_p, -0.5*manif(i).C_p, '-m');
    hold on;
    grid on;
    xlabel('time'); ylabel('Energy E');
    subplot(5,1,4)
    plot(manif(i).time_s_p, manif(i).rr_manif_norm_s_p, '-m');
    hold on;
    grid on;
    xlabel('time'); ylabel('r from m1');
    subplot(5,1,5)
    plot(manif(i).time_s_p, manif(i).inc_arr, '-m');
    hold on;
    grid on;
    xlabel('time'); ylabel('inclination [deg]');
   
    % append to XY plot of manifolds
    figure(101)
    hold on
    % plot manifolds
    plot(manif(i).rr_s_p(:,1),manif(i).rr_s_p(:,2),'-m');
    plot(manif(i).rr_s_m(:,1),manif(i).rr_s_m(:,2),'-g');
    % plot cloest approach of manifolds
    plot(manif(i).rr_s_p(manif(i).rr_minIndx_s_p,1),manif(i).rr_s_p(manif(i).rr_minIndx_s_p,2),'xb');
    grid on;
    xlabel('x'); ylabel('y');
    
    % 3D plot of manifolds
    figure(102)
    % plot manifolds
    plot3(manif(i).rr_s_p(:,1),manif(i).rr_s_p(:,2),manif(i).rr_s_p(:,3),'-m');
    plot3(manif(i).rr_s_m(:,1),manif(i).rr_s_m(:,2),manif(i).rr_s_m(:,3),'-g');
    hold on
    % plot cloest approach of manifolds
    plot3(manif(i).rr_s_p(manif(i).rr_minIndx_s_p,1),manif(i).rr_s_p(manif(i).rr_minIndx_s_p,2),...
        manif(i).rr_s_p(manif(i).rr_minIndx_s_p,3),'xb');
    hold on
    grid on;
    xlabel('x'); ylabel('y');
    
    % append location on 3D plot of Halo
    figure(1)
    hold on;
    plot3(manif(i).X0(1),manif(i).X0(2),manif(i).X0(3),'xb');
    
end

% over-plot GEO altitude for reference
figure(101)
hold on;
plotCircle(0,0,(35786+6378)/Lstar, '--', 'k');


%% Find closest approach among all manifolds (in terms of rr)
for i = 1:length(manif)
    manif_s_p_min_array(i) = manif(i).rr_min_s_p;
end
[best_rrManif_s_p_rr, best_rrManif_s_p_indx] = min(manif_s_p_min_array);
fprintf('Manifold with closest approach: %5.0f\n',best_rrManif_s_p_indx);

% append best manifold (in terms of rr)
% Halo plot
figure(1)
hold on;
plot3(manif(best_rrManif_s_p_indx).X0(1),manif(best_rrManif_s_p_indx).X0(2),...
    manif(best_rrManif_s_p_indx).X0(3),'xc','LineWidth',1.5);  
% XY plot of manifolds
figure(101)
hold on
plot(manif(best_rrManif_s_p_indx).rr_s_p(:,1),manif(best_rrManif_s_p_indx).rr_s_p(:,2),'-c');
plot(manif(best_rrManif_s_p_indx).X0(1),manif(best_rrManif_s_p_indx).X0(2),...
    'xc','MarkerSize',6,'LineWidth',1.5); 
plot(manif(best_rrManif_s_p_indx).rr_s_p(manif(best_rrManif_s_p_indx).rr_minIndx_s_p,1),...
     manif(best_rrManif_s_p_indx).rr_s_p(manif(best_rrManif_s_p_indx).rr_minIndx_s_p,2),...
    '^c','MarkerSize',6,'LineWidth',1.5); 

% 3D plot of manifolds
figure(102)
hold on
plot3(manif(best_rrManif_s_p_indx).rr_s_p(:,1),manif(best_rrManif_s_p_indx).rr_s_p(:,2),manif(best_rrManif_s_p_indx).rr_s_p(:,3),'-c');
hold on
plot3(manif(best_rrManif_s_p_indx).X0(1),manif(best_rrManif_s_p_indx).X0(2),...
    manif(best_rrManif_s_p_indx).X0(3),'xc','MarkerSize',6,'LineWidth',1.5); 
hold on
plot3(manif(best_rrManif_s_p_indx).rr_s_p(manif(best_rrManif_s_p_indx).rr_minIndx_s_p,1),...
      manif(best_rrManif_s_p_indx).rr_s_p(manif(best_rrManif_s_p_indx).rr_minIndx_s_p,2),...
      manif(best_rrManif_s_p_indx).rr_s_p(manif(best_rrManif_s_p_indx).rr_minIndx_s_p,3),...
      '^c','MarkerSize',6,'LineWidth',1.5); 

% print result location state for patching with manifold
fprintf('Manifold patch point:\n   %1.8f  %1.8f  %1.8f\n   %1.8f  %1.8f  %1.8f\n',...
    manif(best_rrManif_s_p_indx).rr_s_p(manif(best_rrManif_s_p_indx).rr_minIndx_s_p,:),...
    manif(best_rrManif_s_p_indx).vv_s_p(manif(best_rrManif_s_p_indx).rr_minIndx_s_p,:));

    
%% find corresponding location on the halo
bestpatch_indx = manif(best_rrManif_s_p_indx).X0index;
figure(2)
xline(time_halo1(bestpatch_indx)/Thalo1,'--r');
figure(3)
subplot(3,1,1)
xline(time_halo1(bestpatch_indx)/Thalo1,'--r');
subplot(3,1,2)
xline(time_halo1(bestpatch_indx)/Thalo1,'--r');
subplot(3,1,3)
xline(time_halo1(bestpatch_indx)/Thalo1,'--r');








