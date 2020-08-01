function [] = plot_evol_orbit(time, dynmat)
% function plots orbit based on dynamics iterated for orbit raising problem

% extract states
r = dynmat(:,1);
theta = dynmat(:,2);
vr = dynmat(:,3);
vtheta = dynmat(:,4);

% state evolution plot
figure(501)
clf    % clear figure
subplot(2,1,1)
hold on
plot(time,r,'-xb');
plot(time,theta,'-xr');
plot(time,vr,'-xc');
plot(time,vtheta,'-xm');
legend('r','\theta','v_r','v_{\theta}','location','northwest')
xlabel('time');
grid on; grid minor;


% orbit plot
nsteps = length(time);
for i = 1:nsteps
    xy_traj(i,:) = r(i,1) * [cos(theta(i,1)), sin(theta(i,1))];
end
% initial and final positions
theta_ref = linspace(0,2*pi,nsteps)';
for i = 1:nsteps
    xy_initial(i,:) = r(1,1) * [cos(theta_ref(i,1)), sin(theta_ref(i,1))];
    xy_final(i,:) = r(nsteps,1) * [cos(theta_ref(i,1)), sin(theta_ref(i,1))];
end
subplot(2,1,2)
hold on
% plot transfer trajectory
plot(xy_traj(:,1), xy_traj(:,2),'-x');
plot(xy_traj(1,1), xy_traj(1,2),'*');  % marker at initial position
plot(xy_traj(nsteps,1), xy_traj(nsteps,2),'^'); % marker at final position
% initial and final trajectory
plot(xy_initial(:,1), xy_initial(:,2),'-.k');
plot(xy_final(:,1), xy_final(:,2),'--k');
% control vector
% quiver(xy_traj(:,1), xy_traj(:,2),-ux, -uy)
caption_finalorbit = strcat('@r = ', num2str(r(nsteps,1)));
legend('transfer','start','end','@r = 1',caption_finalorbit,...
    'location','southwest')
axis equal
grid on; grid minor;
xlabel('x'); ylabel('y');


end
