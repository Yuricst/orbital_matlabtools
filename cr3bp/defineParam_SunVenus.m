%% Define CR3BP for Sun-Venus parameters
% Yuri Shimane, 2020/04/09
fprintf('Sun-Venus system \n');

systemname = 'SunVenus';
% define CR3BP parameters
m1_dim = 1.98847e30; % Earth mass [kg]
m2_dim = 4.8675e24;  % Moon mass [kg]
mu = m2_dim/(m1_dim+m2_dim);
mstar = m1_dim + m2_dim;
Lstar = 108208000;  % [km] or 384400
G = 6.674*10^-11;      % [m^3/kg.s^2]
Tstar = sqrt((Lstar*1000)^3/(G*mstar));
LP = lagrangePoints(mu); % Lagrange points

fprintf('System mu: %2.6s\n',mu)
fprintf('Characteristic length scale: %6.5s [km]\n',Lstar);
fprintf('Characteristic time scale: %f [days]\n',Tstar/(60*60*24));


