function [V] = get_potential_energy(r1,r2,r3,r4,theta2,theta20,theta30,K2,K3)

% Section for theta2 vector
delta = sqrt(r1^2 + r2^2 - 2.*r1.*r2.*cos(theta2));
beta = acos((r1^2 + delta.^2 - r2^2)./(2.*r1.*delta));
psi = acos((r3^2 + delta.^2 - r4^2)./(2.*r3.*delta));
lamda = acos((r4^2 + delta.^2 - r3^2)./(2.*r4.*delta));
theta3 = psi - beta;
theta4 = pi - lamda - beta;

theta40 = theta4(1);



% psi1 = theta2 - theta20;
psi2 = theta2 - theta20 - (theta3 - theta30);
psi3 = theta4 - theta40 - (theta3 - theta30);

V = 0.5*(K2*psi2.^2 + K3*psi3.^2);




end