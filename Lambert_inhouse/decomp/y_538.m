function y = y_538(r1,r2,A,z)
% equation 538 from Curtis textbook
% input: y_538(r1,r2,A,z)

y = norm(r1) + norm(r2) + A*(z*stumpff_S(z) - 1)/sqrt(stumpff_C(z));

end