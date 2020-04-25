%% Define CR3BP for Sun-Earth parameters
% Yuri Shimane, 2020/03/08
fprintf('Sun-Earth system \n');

systemname = 'SunEarth';
% define CR3BP parameters
mu = 3.040357143e-6;
Lstar = 1.495978714*10^8;  % [km]
G = 6.674*10^-11;      % [m^3/kg.s^2]
mstar = 1.98847e30;    % Sun mass + Earth mass... [kg]
Tstar = sqrt((Lstar*1000)^3/(G*mstar));
LP = lagrangePoints(mu); % Lagrange points

fprintf('System mu: %2.6s\n',mu)
fprintf('Characteristic length scale: %6.5s [km]\n',Lstar);
fprintf('Characteristic time scale: %f [days]\n',Tstar/(60*60*24));



