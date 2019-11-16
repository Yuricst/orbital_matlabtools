%% Keplerian orbital elements
% Display keplerian orbital elements to study 3D orbits
% Yuri Shimane, 2019.11.16
% ========================================= %

% house keeping
clear; close all; clc;

% ======= INPUT orbital elements ======= %
% semi-major axis [km]
a = 13000;
% eccentricity
e = 0.4;
% inclination [deg]
inc = 42;
% RAAN [deg]
RAAN = 65;
% argument of periapsis [deg]
omega = 60;
% true anomaly [deg]
theta = 72;
% center body radius [km]
R_center = 6378; % Earth radius
% ====================================== %
if e > 1
    disp('Eccentricity should be less than 1 for closed trajectory')
    return
end

% gravitational parameter
mu = 398600; % [km^3/s^2]

% orbital period
P = 2*pi()*sqrt(a^3/mu);

% periapsis distance
rp = a*(1-e);

% convert elements to state-vector
[r,v] = el2sv(a,e,inc,RAAN,omega,theta,mu);
% go back to orbital elements, where struct OE includes h and e in vector
OE = sv2el(r,v,mu);
% define periapsis vector
rp_vect = (rp/norm(OE.e)) * OE.e;

% specific angular momentum scalar
h = norm(OE.h);

% obtain orbit
n = 200;  % step-size for plotting
theta_range = linspace(0,360,n)'; % range of true-anomaly for plotting
% initialize
r_PF  = zeros(3,n);
r_GEC = zeros(3,n);
for i = 1:n
    % vector in perifocal frame
    r_PF(:,i) = h^2/(mu*(1 + e*cosd(theta_range(i,1))))...
        *[cosd(theta_range(i,1)); sind(theta_range(i,1)); 0];
    % vector in GEC
    tmp1 = RotMat3(-omega,0) * r_PF(:,i);
    tmp2 = RotMat1(-inc,0) * tmp1;
    r_GEC(:,i) = RotMat3(-RAAN,0) * tmp2;
end

% obtain spacecraft location
r_PF_sc = h^2/(mu*(1 + e*cosd(theta)))...
    *[cosd(theta); sind(theta); 0];
tmp1 = RotMat3(-omega,0) * r_PF_sc;
tmp2 = RotMat1(-inc,0) * tmp1;
r_GEC_sc(:,1) = RotMat3(-RAAN,0) * tmp2;


%% Plot results
% plotting parameters
lw = 0.7;
fs = 12;

% create 2D-object
[x_circle,y_circle] = circleplot(0,0,R_center);

% PF-frame plot (2D)
figure(11)
plot(r_PF(1,:),r_PF(2,:));
hold on
plot(r_PF_sc(1,1),r_PF_sc(2,1),'*','markersize',5)
hold on
plot(x_circle,y_circle,'b')
title('Orbit in perifocal frame')
axis equal
grid on; grid minor;
xline(0,'--k','linewidth',lw);
yline(0,'--k','linewidth',lw);
xlabel('x_{PF} [km]');
ylabel('y_{PF} [km]');
set(gca,'fontsize',fs);

% create 3D-center object


% GEC-frame plot (3D)
figure(22)
plot3(r_GEC(1,:),r_GEC(2,:),r_GEC(3,:))
hold on
plot3(r_GEC_sc(1,1),r_GEC_sc(2,1),r_GEC_sc(3,1),'*r','markersize',5)
hold on
quiver3(0,0,0,...
    r_GEC_sc(1,1),r_GEC_sc(2,1),r_GEC_sc(3,1),'--r')
hold on
plot3(0,0,0,'.b','markersize',128)
hold on
quiver3(0,0,0,...
    rp_vect(1,1), rp_vect(2,1), rp_vect(3,1),'--') % periapsis vector
quiver3(0,0,0,...
    0.2*OE.h(1,1), 0.2*OE.h(2,1), 0.2*OE.h(3,1)) % specific h vector
% vernal equinox
quiver3(0,0,0,...
    1.75*R_center,0,0,'--k','linewidth',1)
% RAAN vector
quiver3(0,0,0,...
    1.75*R_center*cosd(RAAN),1.75*R_center*sind(RAAN),0,'--k','linewidth',1)
% orbital zenith
quiver3(0,0,0,...
    0,0,1.75*R_center,'--k','linewidth',1)
xlabel('x_{GEC} [km]');
ylabel('y_{GEC} [km]');
zlabel('z_{GEC} [km]');
title('Orbit in GEC')
axis equal
grid on; grid minor;
set(gca,'fontsize',fs);





