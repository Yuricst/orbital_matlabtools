function S = stumpff_S(z)
% Stumpff S(z) function
if z > 0
    S = ( sqrt(z) - sin(sqrt(z)) ) / z^1.5;
elseif z == 0
    S = 1/6;
elseif z < 0
    S = ( sinh(sqrt(-z) - sqrt(-z)) ) / (-z)^1.5;
end
end
