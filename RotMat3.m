function R3 = RotMat3(phi,unit)
% Rotatonal matrix about third axis
% Input: RotMat1(phi,unit)
%   phi (angle)
%   unit : deg (0) or radian (1)
% Output: R3
%   3 by 3 rotational matrix

if unit == 0
    %compute matrix in degree
    R3 = [cosd(phi) sind(phi) 0;...
        -sind(phi) cosd(phi) 0;...
        0 0 1];
else %(if unit == 1)
    %comute matrix in radian
    R3 = [cos(phi) sin(phi) 0;...
        -sin(phi) cos(phi) 0;...
        0 0 1];
end

end