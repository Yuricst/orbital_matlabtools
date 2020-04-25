function [v0,vf] = lambert(r0,rf,tof,mu,grade)
% Curtis formulation of Lambert problem (TPBVP) 
% Matlab transcription based on textbook-mfiles but debugged
% INPUT
%   r0 : initial position 1x3 array [km]
%   rf : final position 1x3 array [km]
%   tof : time of flight [sec]
%   mu : gravitational parameter [km^3/s^2]
%   grade : (string) 'prograde' or 'retrograde'
% OUTPUT
%   v0 : initial velocity 1x3 array [km/s]
%   vf : final velocity 1x3 array [km/s]
% ====================================================== %

% compute angular separation
dtheta = acos(dot(r0,rf)/(norm(r0)*norm(rf)));
c12 = cross(r0,rf);

% update dtheta for special cases
if strcmp(grade,'retrograde') == 1
    if c12(2) >= 0
        dtheta = 2*pi - dtheta;
    end
else
    if c12(2) <= 0
        dtheta = 2*pi - dtheta;
    end
end

% compue input parameter A
A = sin(dtheta) * sqrt( norm(r0)*norm(rf) / (1-cos(dtheta)) );
fprintf('A = %f\n',A)

% find initial value of z
z = -100;
if fobj(r0,rf,A,z,tof,mu) < 0
    fprintf('z = %f, F = %f\n',z,fobj(r0,rf,A,z,tof,mu))
    z = z + 0.1;
end
fprintf('Initial guess of z: %f\n',z)

% set error tolerance and limit of number of iterations
tol   = 1.e-8;
nmax  = 5000;
% initialize
ratio = 1;
n = 1;
% solve with newton-raphson method
while abs(ratio) > tol
    fprintf('Current iteration: %d\n',n)
    fprintf('F: %f\n', fobj(r0,rf,A,z,tof,mu))
    fprintf('dFdz: %f\n',dfobj(r0,rf,A,z))
    % compute ratio
    ratio = fobj(r0,rf,A,z,tof,mu) / dfobj(r0,rf,A,z);
    % iterate value of z with newton-raphson
    z     = z - ratio;
    fprintf('z: %f\n',z)
    if n > nmax
        fprintf('\n\n **Number of iterations exceeds %g in ''lambert'' \n\n ',nmax)
        break
    end
    n = n + 1;
end

% formulation below is same as Curtis
%...Equation 5.46a:
f    = 1 - y_538(r0,rf,A,z)/norm(r0);
%...Equation 5.46b:
g    = A*sqrt(y_538(r0,rf,A,z)/mu);
%...Equation 5.46d:
gdot = 1 - y_538(r0,rf,A,z)/norm(rf);

%...Equation 5.28:
v0   = 1/g*(rf - f*r0);
%...Equation 5.29:
vf   = 1/g*(gdot*rf - r0);

% --------------- nested functions --------------- %
    function y = y_538(r1,r2,A,z)
        % equation 538 from Curtis textbook
        % input: y_538(r1,r2,A,z)

        y = norm(r1) + norm(r2) + A*(z*stumpff_S(z) - 1)/sqrt(stumpff_C(z));

    end

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

    function fval = fobj(r1,r2,A,z,tof,mu)
        % objective funtion
        % input : fobj(r1,r2,A,z,tof,mu)
        fval = (y_538(r1,r2,A,z)/stumpff_C(z))^1.5*stumpff_S(z) + A*sqrt(y_538(r1,r2,A,z)) - sqrt(mu)*tof;
    end

    function dfval = dfobj(r1,r2,A,z)
        % objective funtion
        % input: dfobj(r1,r2,A,z)
        if z == 0
            dfval = sqrt(2)/40*y_538(r1,r2,A,0)^1.5 + A/8*(sqrt(y_538(r1,r2,A,0)) + A*sqrt(1/2/y_538(r1,r2,A,0)));
        else
            dfval = (y_538(r1,r2,A,z)/stumpff_C(z))^1.5*(1/2/z*(stumpff_C(z) - 3*stumpff_S(z)/2/stumpff_C(z)) ...
                   + 3*stumpff_S(z)^2/4/stumpff_C(z)) + A/8*(3*stumpff_S(z)/stumpff_C(z)*sqrt(y_538(r1,r2,A,z)) ...
                   + A*sqrt(stumpff_C(z)/y_538(r1,r2,A,z)));
        end
    end
% ------------------------------------------------ %
end

        
        