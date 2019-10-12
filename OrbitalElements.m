function OE = OrbitalElements(r,v,mu)
% COMPUTE ORBITAL ELEMENTS from r, v, mu
% Yuri SHIMANE, 2018.10.19
%   position vector r = [rx,ry,rz] (km)
%   velocity vector r = [vx,vy,vz] (km/s)
%   gravitational parameter mu (km^3.s^-2)
% OUTPUT: 
% Structure array
%   OE.h  : specific angular momentum vector (km^2.s^-1)
%   OE.i  : inclination angle (deg)
%   OE.RA : right ascension (deg)
%   OE.e  : eccentricity vector
%   OE.omega : argument of perigee (deg)
%   OE.theta : true anomaly (deg)

%specific angular momentum vector
OE.h = cross(r,v);

%inclination angle (deg)
OE.i = acosd(OE.h(1,3)/norm(OE.h));

%eccentricity vector
OE.e = cross(v,OE.h)/mu - r/norm(r);

%right ascension
K = [0,0,1];    %GEC z-axis
N = cross(K,OE.h);
if N(1,2) > 0   %if Ny > 0 then RA is the cosine
    OE.RA = acosd(N(1,1)/norm(N));   %right ascension
else            %elseif Ny < 0 RA = 360 - cosine
    OE.RA = 360 - acosd(N(1,1)/norm(N));   %right ascension
end

%argument of perigee
if OE.e(1,3) > 0   %if ez > 0 then omega < 180
    OE.omega = acosd(dot(OE.e,N)/(norm(OE.e)*norm(N)));
else            %elseif ez < 0 then omega > 180
    OE.omega = 360 - acosd(dot(OE.e,N)/(norm(OE.e)*norm(N)));
end

%true anomaly
v_r = dot(v,r)/norm(r);  %radial velocity
if v_r > 0  %if v_r > 0, theta<180 (moving away)
    OE.theta = acosd(dot(OE.e,r)/(norm(OE.e)*norm(r)));
else        %if v_r < 0, theta>180 (coming back)
    OE.theta = 360 - acosd(dot(OE.e,r)/(norm(OE.e)*norm(r)));
end

% K = [0,0,1];    %GEC z-axis
% N = cross(K,OE.h);
% OE.RA = acosd(N(1,1)/norm(N));   %right ascension
% OE.omega = acosd(dot(OE.e,N)/(norm(OE.e)*norm(N)));  %argument of perigee
% OE.theta = acosd(dot(OE.e,r)/(norm(OE.e)*norm(r)));  %true anomaly
% 

end
