function out = legendre_cn(n,mu,gammaL,Lpoint)
% function comutes Legendre constant c_n for CR3BP system in colinear
% system
% FORMAT : out = legendre_cn(n,mu,gammaL)
% INPUT :
%   n : >1 (?)
%   mu : CR3BP system mu
%   gammaL : distance to L point of interest
%   Lpoint : 1, 2, or 3
% NOTE : Validated using Richardson (1980)
% Yuri Shimane, 2020/03/05
% ============================================================== 

if Lpoint == 1
    out = (1/gammaL^3) * (mu + (-1)^n * (1-mu)*gammaL^(n+1) / (1 - gammaL)^(n+1));
elseif Lpoint ==2
    out = (1/gammaL^3) * ((-1)^n*mu + (-1)^n * (1-mu)*gammaL^(n+1) / (1 + gammaL)^(n+1));
elseif Lpoint ==3
    out = (1/gammaL^3) * (1 -mu + mu*gammaL^(n+1) / (1 + gammaL)^(n+1));
end

end