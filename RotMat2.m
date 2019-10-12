function R2 = RotMat2(phi,unit)
% Rotatonal matrix about second axis
% Input: RotMat1(phi,unit)
%   phi (angle)
%   unit : deg (0) or radian (1)
% Output: R3
%   3 by 3 rotational matrix

if unit == 0
    %compute matrix in degree
    R2 = [cosd(phi) 0 -sind(phi);...
        0 1 0;...
        sind(phi) 0 cosd(phi)];
else %(if unit == 1)
    %comute matrix in radian
    R2 = [cos(phi) 0 -sin(phi);...
        0 1 0;...
        sin(phi) 0 cos(phi)];
end

end