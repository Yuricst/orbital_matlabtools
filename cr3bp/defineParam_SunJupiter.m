%% Define CR3BP for Earth-Moon parameters
% Yuri Shimane, 2020/05/23
fprintf('Sun-Jupiter system \n');

systemname = 'SunJupiter';
% define CR3BP parameters
mu = 9.537e-4;
%mstar = m1_dim + m2_dim;
Lstar = 7.784e8;  % [km] or 384400
G = 6.674*10^-11;      % [m^3/kg.s^2]
% mstar = 1.98847e30;    % Sun mass + Earth mass... [kg]
Tstar = 3.733e8 / (2*pi);
%Tstar = sqrt((Lstar*1000)^3/(G*mstar));
LP = lagrangePoints(mu); % Lagrange points

fprintf('System mu: %2.6s\n',mu)
fprintf('Characteristic length scale: %6.5s [km]\n',Lstar);
fprintf('Characteristic time scale: %f [days]\n',Tstar/(60*60*24));


