function [lambda_val, lambda_realarr, lambda_imagarr] = linear_cr3bp_lambda(c2)
% funtion computes roots of characteristic equation for coupled motion in
% XY plane
% FORMAT: [lambda, lambda_realarr, lambda_imagarr] = linear_cr3bp_lambda(c2)
% Yuri Shimane, 2020/03/07
% ========================================================== %

% initialize
lambda_val = [];

syms lambda 
f = lambda^4 + (c2-2)*lambda^2 - (c2 - 1)*(1 + 2*c2);
% solve function
S = solve(f,0);
j = 1; k = 1;
for i = 1:4
    if imag(double(S(i,1))) == 0
        lambda_realarr(j,1) = double(S(i,1));
        j = j+1;
        % return value of lambda
        lambda_val = abs(lambda_realarr(1,1));
    else
        lambda_imagarr(k,1) = double(S(i,1));
        k = k+1;
    end
end

end