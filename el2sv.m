function [r,v] = el2sv(a,e,inc,RAAN,omega,theta,mu)
% COMPUTE STATE-VECTOR FROM ORBITAL ELEMENTS
% Yuri SHIMANE, 2019.11.15
% INPUT:
%   a : semi-major axis (km)
%   inc : inclination (deg)
%   RAAN : Right Ascension of Ascending Node (deg)
%   omega : argument of periapsis (deg)
%   e : eccentricity scalar
%   theta : true anomaly (deg)
%   mu : gravitational parameter (km^3.s^-2)
% OUTPUT: 
%   r : 3x1 position vector in GEC (km) 
%   v : 3x1 velocity vector in GEC (km/s)
% ========================================================= %

% compute specific angular momentum
h = sqrt(a*mu*(1-e^2));

% perifocal position vector
rPF = h^2/(mu*(1+e*cosd(theta))) * [cosd(theta); sind(theta); 0];
vPF = (mu/h) * [-sind(theta); e+cosd(theta); 0];

% define rotational matrices for 1st axis
function R1 = RotMat1(phi)
% Rotatonal matrix about first axis
% Input: RotMat1(phi,unit)
%   phi (angle) [deg]
% Output: R3
%   3 by 3 rotational matrix
    R1 = [1 0 0;...
        0 cosd(phi) -sind(phi);...
        0 -sind(phi) cosd(phi)];
end

% define rotational matrices for 3rd axis
function R3 = RotMat3(phi)
% Rotatonal matrix about third axis
% Input: RotMat1(phi,unit)
%   phi (angle) [deg]
% Output: R3
%   3 by 3 rotational matrix
    R3 = [cosd(phi) sind(phi) 0;...
        -sind(phi) cosd(phi) 0;...
        0 0 1];
end

% compute position in GEC
rtmp1 = RotMat3(-omega) * rPF;
rtmp2 = RotMat1(-inc) * rtmp1;
r     = RotMat3(-RAAN) * rtmp2;
% compute velocity in GEC
vtmp1 = RotMat3(-omega) * vPF;
vtmp2 = RotMat1(-inc) * vtmp1;
v     = RotMat3(-RAAN) * vtmp2;

end