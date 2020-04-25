function inc_arr = inclinationArray(hh)
% Function computes array of inclination
% Yuri Shimane, 2020/04/25
% FORMAT
%   inc_arr = inclinationArray(hh)
% INPUT
%   hh : n by 3 array of angular momentum hx, hy, hz
% OUTPUT
%   inc_arr : n by 1 array of inclination in degrees
% ================================================ %

[n,~] = size(hh);
inc_arr = zeros(n,1);

for i = 1:n
    inc_arr(i,1) = acosd(hh(i,3)/norm(hh(i,:)));
end

end