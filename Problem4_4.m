function [theta3,theta4,xp,yp,MA] = Problem4_4(r1,r2,r3,r4,a3,b3)
% Develop a program to calculate the position of a four-bar mechanism

% theta2 is the input
theta2 = 0:10:360;

theta1 = zeros(size(theta2));
theta3 = zeros(size(theta2));
theta4 = zeros(size(theta2));
xp = zeros(size(theta2));
yp = zeros(size(theta2));
MA = zeros(size(theta2));

% Calculate relative crank angle
theta2pr = theta2 - theta1;

delta = (r1^2 + r2^2 - 2*r1*r2*cosd(theta2pr)).^0.5;

beta = acosd((r1^2 + delta.^2 - r2^2)./(2*r1*delta));

psi = acosd((r3^2 + delta.^2 - r4^2)./(2*r3*delta));

lamda = acosd((r4^2 + delta.^2 - r3^2)./(2*r4*delta));

for i = 1:length(theta2pr)
    if theta2pr(i) >= 0 && theta2pr(i) <= 180
        theta3(i) = psi(i) - (beta(i) - theta1(i));
        theta4(i) = 180 - lamda(i) - (beta(i) - theta1(i));
    else
        theta3(i) = psi(i) + (beta(i) + theta1(i));
        theta4(i) = 180 - lamda(i) + (beta(i) + theta1(i));
    end
end

MA = (r4*sind(theta4 - theta3))./(r2*sind(theta2 - theta3));

xp = r2*cosd(theta2) + a3*cosd(theta3) - b3*sind(theta3);
yp = r2*sind(theta2) + a3*sind(theta3) + b3*cosd(theta3);


%% Plot results
% Plot relative angles theta3 and theta4 against theta2, the input angle
figure(1);
plot(theta2,theta3);
hold on
plot(theta2,theta4);
hold off
xlabel('\theta_2 (deg)');
ylabel('Relative angle (deg)');
legend('\theta_3','\theta_4');
title('Relative Angles vs. Input Angle');

% Plot motion path of coupler point
figure(2);
plot(xp,yp);
title('Path of Coupler Point P');

% Plot mechanical advantage against input angle
figure(3);
plot(theta2,MA);
xlabel('\theta_2 (deg)');
ylabel('Mechanical Advantage');
title('Mechanical Advantage vs. Input Angle');




end