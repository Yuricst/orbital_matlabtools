%% Define CR3BP for Kathleen Howell (1984) Table 1
% Yuri Shimane, 2020/03/08
fprintf('Howell (1984) Table 1\n');

systemname = 'Howell1984';
% define CR3BP parameters
mu = 0.04;
mstar = 1;
Lstar = 1;  % [km] or 384400
G = 6.674*10^-11;      % [m^3/kg.s^2]
% mstar = 1.98847e30;    % Sun mass + Earth mass... [kg]
Tstar = sqrt((Lstar*1000)^3/(G*mstar));
LP = lagrangePoints(mu); % Lagrange points

fprintf('System mu: %2.6s\n',mu)
fprintf('Characteristic length scale: %6.5s [km]\n',Lstar);
fprintf('Characteristic time scale: %f [days]\n',Tstar/(60*60*24));


% Table 1
Howell84_table1 = [...
    0.723268 0.729988 0.753700 0.777413 0.801125 0.817724 ;
    0.040000 0.215589 0.267595 0.284268 0.299382 0.313788 ;
    0.198019 0.397259 0.399909 0.361870 0.312474 0.271306 ;
    1.300177 1.348532 1.211253 1.101099 1.017241 0.978635 ];

[~, tablew] = size(Howell84_table1);

