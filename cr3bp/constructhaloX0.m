function [X0_halo,Thalo,stm_T2,X0_analytical] = constructhaloX0(mu,Lstar,phi,lp,northsouth,Az_km)
% Function to create Halo orbit by grouping other functions
% Yuri Shimane, 2020/04/17
% =============================================== %

% analytical 3rd order solution
[X0_analytical,T_analytical,~,~,~] = ...
    collinearHalo_analytical(mu,Lstar,northsouth,lp,Az_km,phi);

% Differential-correction using event condition in ode45
% tolerance of residuals for DC process
residualtol = 1e-10;
% rough estimate for integration number of steps
nsteps = 5000;
% tolerance setting for ode45
reltol = 1e-10;
abstol = 1e-10;
% run DC algorithm
[X0_halo,Thalo,~,~,~,stm_T2,~,~] = ...
    collinearHalo_differentialCorrection(mu,X0_analytical,T_analytical/2,...
                                            residualtol,nsteps,reltol,abstol);

end