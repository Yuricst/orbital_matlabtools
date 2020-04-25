function [] = plotCircle(centerX, centerY, r, pform, pcol)
% Function plots a circle centered about centerX and centerY with radius r

% coordinates
theta = linspace(0, 2*pi, 500);
x = zeros(1,length(theta));
y = zeros(1,length(theta));
for i = 1:length(theta)
    x(i) = r * cos(theta(i)) + centerX;
    y(i) = r * sin(theta(i)) + centerY;
end
plot(x,y,pform,'Color',pcol);

end