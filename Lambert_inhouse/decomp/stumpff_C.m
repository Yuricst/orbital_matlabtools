function C = stumpff_C(z)
% Stumpff C(z) function
if z > 0
    C = ( 1 - cos(sqrt(z)) ) / z;
elseif z == 0
    C = 1/2;
elseif z < 0
    C = (cosh(sqrt(-z)) - 1)/(-z);
end
end