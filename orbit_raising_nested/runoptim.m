function [u,fval,exitflag,output,cLast,ceqLast] = runoptim(u0,opts,nsteps,mcxop)
% function calling fmincon, with embedded objective and nlncon functions
% =================================================== %
% Yuri Shimane, 2020/02/20
% INPUT
%   u0 : initial guess, concatenated (nsteps x 1) vector, [u1 u2] ([ur utheta])
%   opts : optimoptions handle
%   nsteps : number of steps to discretize trajectory (e.g. 100)
%   mcxop : option to calculate gradient with MCX ('true' / 'false')
% =================================================== %

if nargin == 1 % no options supplied
    opts = [];
end

% declare global variables for preventing unecessary dynamics computation
uLast = [];
fvalLast = [];
gradfLast = [];
cLast = [];
ceqLast = [];
gradcLast = [];
gradceqLast = [];

% final time
tf = 3.32;

% objective function, nested below
ofun = @objfun;
% constraint function, nested below
cfun = @nlcon;
% lower and upper bounds
lb = [-ones(nsteps,1);
      -ones(nsteps,1)];
ub = [ones(nsteps,1);
      ones(nsteps,1)];

% call fmincon
[u,fval,exitflag,output] = fmincon(ofun,u0,[],[],[],[],lb,ub,cfun,opts);

% ... nested functions below ...
    % ===== objective function ===== %
    function [fval,gradf] = objfun(u)
%         disp('running objective function...')
        % check if computation is necessary
        if ~isequal(u,uLast)
            % if u is not equal to uLast, run dynamics
            [fvalLast,gradfLast,cLast,ceqLast,gradcLast,gradceqLast] = rundynamics(u,nsteps,tf,mcxop);
            uLast = u;
        end
        % return objective function value and store updated value
        fval = fvalLast;
        gradf = gradfLast;   
    end

    % ===== non-linear constraint function ===== %
    function [c,ceq,gradc,gradceq] = nlcon(u)
%         disp('running nlcon...')
        % check if computation is necessary
        if ~isequal(u,uLast)
            % if u is not equal to uLast, run dynamics
            [fvalLast,gradfLast,cLast,ceqLast,gradcLast,gradceqLast] = rundynamics(u,nsteps,tf,mcxop);
            uLast = u;
        end
        % return objective function value and store updated value
        c = cLast;
        ceq = ceqLast;
        gradc = gradcLast;
        gradceq = gradceqLast;
    end
    
end


