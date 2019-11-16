function [x,y] = circleplot(x0,y0,r)
% plot circle
theta = linspace(0,360,1000);

x = x0 + r*cosd(theta);
y = y0+ r*sind(theta);

end

