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