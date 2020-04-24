function [a21,a22,a23,a24,a31,a32,b21,b22,b31,b32,d21,d31,d32,...
          d1,d2,a1,a2,s1,s2,l1,l2] = legendre_coefficients(c2,c3,c4,lambda,k)
% funtion computes coefficients for third-order solution in cr3bp
% FORMAT: out = legendre_coefficients(c2,lambda,k)
% Yuri Shimane, 2020/03/07
% ========================================================== %

% coefficient d1, d2
d1 = ( 3*(lambda^2)/k ) * (k*(6*lambda^2 - 1) - 2*lambda);
d2 = ( 8*(lambda^2)/k ) * (k*(11*lambda^2 - 1) - 2*lambda);
% coefficient d21
d21 = -c3/(2*lambda^2);
% coefficient b21, b22
b21 = ( -3*c3*lambda/(2*d1) )*(3*k*lambda - 4);
b22 = 3*c3*lambda/d1;
% coefficient a21 ~ a24
a21 = 3*c3*(k^2 - 2)/(4*(1+2*c2));
a22 = 3*c3/(4*(1+2*c2));
a23 = (-3*c3*lambda/(4*k*d1)) * (3*k^3*lambda - 6*k*(k - lambda) + 4);
a24 = (-3*c3*lambda/(4*k*d1)) * (2 + 3*k*lambda);
% b31 ~ b32
b31 = (3/(8*d2))*(8*lambda*(3*c3*(k*b21-2*a23)-c4*(2+3*k^2)) ...
      + (9*lambda^2 + 1 + 2*c2)*(4*c3*(k*a23-b21)+k*c4*(4+k^2)));
b32 = (1/d2)*(9*lambda*(c3*(k*b22+d21-2*a24)-c4) + (3/8)*(9*lambda^2 + 1 ...
      + 2*c2)*(4*c3*(k*a24 - b22) + k*c4));
% a31 ~ a32
a31 = (-9*lambda/(4*d2)) * (4*c3*(k*a23 - b21) + k*c4*(4+k^2)) ...
      + ((9*lambda^2 + 1 - c2)/(2*d2))*(3*c3*(2*a23 - k*b21) + c4*(2+3*k^2));
a32 = (-1/d2)*((9*lambda/4)*(4*c3*(k*a24 - b22) + k*c4) ...
      + (3/2)*(9*lambda^2 + 1 - c2)*(c3*(k*b22 + d21 - 2*a24) - c4));
% coefficient d31 ~ d32
d31 = (3/(64*lambda^2))*(4*c3*a24 + c4);
d32 = (3/(64*lambda^2))*(4*c3*(a23 - d21) + c4*(4+k^2));
% frequency correction s1, s2
s1 = (1/(2*lambda*(lambda*(1+k^2)-2*k)))*((3/2)*c3*(2*a21*(k^2-2) ...
     - a23*(k^2+2)-2*k*b21) - (3/8)*c4*(3*k^4-8*k^2+8));
s2 = (1/(2*lambda*(lambda*(1+k^2)-2*k)))*((3/2)*c3*(2*a22*(k^2-2) ...
     + a24*(k^2+2)+2*k*b22+5*d21) + (3/8)*c4*(12-k^2));
% amplitude constraint param. a1, a2, l1, l2
a1 = -(3/2)*c3*(2*a21 + a23 + 5*d21) - (3/8)*c4*(12 - k^2);
a2 = (3/2)*c3*(a24 - 2*a22) + (9/8)*c4;
l1 = a1 + 2*lambda^2 * s1;
l2 = a2 + 2*lambda^2 * s2;

end

