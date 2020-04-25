function fval = fobj(r1,r2,A,z,tof,mu)
% objective funtion
% input : fobj(r1,r2,A,z,tof,mu)

fval = (y_538(r1,r2,A,z)/stumpff_C(z))^1.5*stumpff_S(z) + A*sqrt(y_538(r1,r2,A,z)) - sqrt(mu)*tof;

end