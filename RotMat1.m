function R1 = RotMat1(phi,unit)
% Rotatonal matrix about first axis
% Input: RotMat1(phi,unit)
%   phi (angle)
%   unit : deg (0) or radian (1)
% Output: R3
%   3 by 3 rotational matrix

if unit == 0
    %compute matrix in degree
    R1 = [1 0 0;...
        0 cosd(phi) sind(phi);...
        0 -sind(phi) cosd(phi)];
else %(if unit == 1)
    %comute matrix in radian
    R1 = [1 0 0;...
        0 cos(phi) -sin(phi);...
        0 -sin(phi) cos(phi)];
end

end