%% Define CR3BP for Earth-Moon parameters
% Yuri Shimane, 2020/03/08
fprintf('Earth-Moon system \n');

systemname = 'EarthMoon';
% define CR3BP parameters
m1_dim = 5.9722e24; % Earth mass [kg]
m2_dim = 7.347e22;  % Moon mass [kg]
% mu = m2_dim/(m1_dim+m2_dim);
mu = 0.0121409319;
mstar = m1_dim + m2_dim;
Lstar = 385692.5;  % [km] or 384400
G = 6.674*10^-11;      % [m^3/kg.s^2]
% mstar = 1.98847e30;    % Sun mass + Earth mass... [kg]
Tstar = sqrt((Lstar*1000)^3/(G*mstar));
LP = lagrangePoints(mu); % Lagrange points

fprintf('System mu: %2.6s\n',mu)
fprintf('Characteristic length scale: %6.5s [km]\n',Lstar);
fprintf('Characteristic time scale: %f [days]\n',Tstar/(60*60*24));

% Earth SOI
SOI_Earth = 925000/Lstar;