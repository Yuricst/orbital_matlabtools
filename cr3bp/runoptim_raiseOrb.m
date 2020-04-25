function [u,fval,exitflag,output] = ...
    runoptim_raiseOrb(mu,u0,mass0,X0,Xf,Isp,Fmax,Tfmax,opts,nsteps,Lstar,Tstar)
% function calling fmincon, with embedded objective and nlncon functions
% Yuri Shimane, 2020/04/17
% =================================================== %


% initialize variables for preventing unecessary dynamics computation
uLast = [];
fvalLast = [];
cLast = [];
ceqLast = [];

% lower and upper bounds
lb = [-ones(nsteps,1);
      -ones(nsteps,1);
      -ones(nsteps,1);
      0];
ub = [ones(nsteps,1);
      ones(nsteps,1);
      ones(nsteps,1);
      Tfmax];

% call fmincon
[u,fval,exitflag,output] = fmincon(@objfun,u0,[],[],[],[],lb,ub,@nlcon,opts);

% ... nested functions below ...
    % ===== objective function ===== %
    function fval = objfun(u)
        % check if computation is necessary
        if ~isequal(u,uLast)
            % if u is not equal to uLast, run dynamics
            [fvalLast,cLast,ceqLast] = ...
                rundynamics_raiseOrb(mu,u,mass0,X0,Xf,Isp,Fmax,nsteps,Lstar,Tstar);
            uLast = u;
        end
        % return objective function value and store updated value
        fval = fvalLast; 
    end

    % ===== non-linear constraint function ===== %
    function [c,ceq] = nlcon(u)
        % check if computation is necessary
        if ~isequal(u,uLast)
            % if u is not equal to uLast, run dynamics
            [fvalLast,cLast,ceqLast] = ...
                rundynamics_raiseOrb(mu,u,mass0,X0,Xf,Isp,Fmax,nsteps,Lstar,Tstar);
            uLast = u;
        end
        % return constraints and store updated value
        c = cLast;
        ceq = ceqLast;
    end

end


