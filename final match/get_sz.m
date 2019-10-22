function [ new_sz ] = get_sz(target_sz, theta)
% Ğı×ªºóµÄ³ß´ç

if theta > 90
    theta = 180-theta;
elseif theta < -90
    theta = -180-theta;
end

theta = abs(theta);
theta1 = atand(target_sz(1)/target_sz(2));
theta2 = theta1 + theta;
theta3 = theta - theta1;
r = sqrt((target_sz(1)/2).^2 + (target_sz(2)/2).^2);
h = r * sind(theta2);
w = r * cosd(theta3);
new_sz(1) = 2*h;
new_sz(2) = 2*w;
end

