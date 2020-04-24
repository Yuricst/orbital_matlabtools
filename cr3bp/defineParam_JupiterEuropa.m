%% Define CR3BP for Jupiter-Europa System
% Yuri Shimane, 2020/03/08
fprintf('Jupiter-Europa system \n');

systemname = 'JupiterEuropa';
% define CR3BP parameters
mu = 2.52800e-5;
Lstar = 670900;  % [km]
G = 6.674*10^-11;      % [m^3/kg.s^2]
% mstar = 1.98847e30;    % Sun mass + Earth mass... [kg]
Tstar = 1.769322 * (24*60*60); % [sec]
LP = lagrangePoints(mu); % Lagrange points

fprintf('System mu: %2.6s\n',mu)
fprintf('Characteristic length scale: %6.5s [km]\n',Lstar);
fprintf('Characteristic time scale: %f [days]\n',Tstar/(60*60*24));



