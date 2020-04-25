function [C,Ubar] = jacobiConst(mu,rr,vv)
% Function computes Jacobi constant C = -2E in CR3BP
% INPUT
%   mu : CR3BP mu
%   rr : array of positions (N by 3)
%   vv : array of velocities (N by 3)
% OUTPUT
%   C  : array of Jacobi constant (N by 1)
%   Ubar : array of pseudo-potential (N by 1)
% ========================================= %

% compute size of array (number of time stpes)
[N,~] = size(rr);
C = zeros(N,1);
% mu of each body
mu1 = 1-mu;
mu2 = mu; 

for i = 1:N
    % compute norm of radii from body 1 and from body 2
    r1_vec = [rr(i,1)+mu, rr(i,2), rr(i,3)];
    r2_vec = [rr(i,1)+(1-mu), rr(i,2), rr(i,3)];
    r1 = norm(r1_vec);
    r2 = norm(r2_vec);
    % compute pseudo potential
    Ubar(i,1) = -(1/2)*(mu1*r1^2 + mu2*r2^2) - mu1/r1 - mu2/r2; 
    % extract velocities
    xdot = vv(i,1);
    ydot = vv(i,2);
    zdot = vv(i,3);
    % Jacobi integral
    C(i,1) = -(xdot^2 + ydot^2 + zdot^2) - 2*Ubar(i,1);
    
end

end