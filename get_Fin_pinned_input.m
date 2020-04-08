function [Fin] = get_Fin_pinned_input(r1,r2,r3,r4,r5,theta2,dy,a3,b3,theta20,theta30,K2,K3)

theta1 = 0;
delta = sqrt(r1^2 + r2^2 - 2.*r1.*r2.*cos(theta2));
beta = acos((r1^2 + delta.^2 - r2^2)./(2.*r1.*delta));
psi = acos((r3^2 + delta.^2 - r4^2)./(2.*r3.*delta));
lamda = acos((r4^2 + delta.^2 - r3^2)./(2.*r4.*delta));
theta3 = psi - beta;
theta4 = pi - lamda - beta;

h32 = (r2.*sin(theta4-theta2))./(r3.*sin(theta3-theta4));
h42 = (r2.*sin(theta3-theta2))./(r4.*sin(theta3-theta4));
theta5 = asin((r2.*cos(theta2) + a3.*cos(theta3) - b3.*sin(theta3) + dy)/r5);

T2 = K2.*(theta2 - theta20);
T3 = K3.*(theta3 - theta30);
% 
% Fin = (h32.*(T2+T3) - T2 - T3.*h42)./(r2.*sin(theta5).*cos(theta2) - r2.*cos(theta5).*sin(theta2)...
%     -a3.*h32.*cos(theta5).*sin(theta3) - b3.*h32.*sin(theta5).*sin(theta3)...
%     - b3.*h32.*cos(theta5).*cos(theta3) - a3.*sin(theta5).*cos(theta5));

M2 = 0; % This will need to be a function dependent upon theta2 in the final version

Fin = (h32.*(T2+T3) - h42.*T3 - M2 - T2)./(r2.*sin(theta2) + h32.*a3.*sin(theta3) + h32.*b3.*cos(theta2));


end