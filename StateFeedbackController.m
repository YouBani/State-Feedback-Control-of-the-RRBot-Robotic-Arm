clear; clc; close all;
syms m1 m2 theta1 theta2 r1 r2 l1 l2 I1 I2 theta1_dot theta1_ddot theta2_dot theta2_ddot u1 u2 g 'real'

m1=1; m2=1; l1=1; l2=1 ;r1=0.45; r2=0.45; g=9.81 ;I2= 0.084; I1= 0.084;

% Question a.
eq1= theta1_ddot*(m2*l1^2 + 2*m2*cos(theta2)*l1*r2 + m1*r1^2 + m2*r2^2 + I1 + I2) - theta2_dot*(l1*m2*r2*theta1_dot*sin(theta2) + l1*m2*r2*sin(theta2)*(theta1_dot + theta2_dot)) - u1 + theta2_ddot*(m2*r2^2 + l1*m2*cos(theta2)*r2 + I2) - g*m2*r2*sin(theta1 + theta2) - g*l1*m1*sin(theta1) - g*m1*r1*sin(theta1);
eq2= theta2_ddot*(m2*r2^2 + I2) - u2 + theta1_ddot*(m2*r2^2 + l1*m2*cos(theta2)*r2 + I2) - g*m2*r2*sin(theta1 + theta2) + l1*m2*r2*theta1_dot*sin(theta2)*(theta1_dot + theta2_dot) - l1*m2*r2*theta1_dot*theta2_dot*sin(theta2);

eq1 = subs(eq1,[theta1_dot,theta2_dot,theta1_ddot,theta2_ddot,u1,u2],[0,0,0,0,0,0]);
eq2 = subs(eq2,[theta1_dot,theta2_dot,theta1_ddot,theta2_ddot,u1,u2],[0,0,0,0,0,0]);

sol=solve([eq1==0, eq2==0] ,[theta1,theta2]);

display(sol.theta1)
display(sol.theta2)
% We have 4 equilibrium points:
% [0,0,0,0]
% [pi,0,0,0]
% [0,pi,0,0]
% [pi,pi,0,0]

% Question b.
u=[u1;u2];
x=[theta1,theta2,theta1_dot,theta2_dot];
x1_dot = theta1_dot;
x2_dot= theta2_dot;
x3_dot= (I2*u1 - I2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) + I2*g*l1*m1*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*l1*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
x4_dot= -(I2*u1 - I1*u2 - I2*u2 - l1^2*m2*u2 - m1*r1^2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1^3*m2^2*r2*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + I2*g*l1*m1*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*u1*cos(theta2) - 2*l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*l1*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + g*l1^2*m1*m2*r2*cos(theta2)*sin(theta1) + l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);

X_dot = [x1_dot;x2_dot;x3_dot;x4_dot];

A = jacobian(X_dot, x);
B = jacobian(X_dot, u);

A1 = subs(A,[theta1,theta2,theta1_dot,theta2_dot],[0,0,0,0]);
B1 = subs(B,[theta1,theta2,theta1_dot,theta2_dot],[0,0,0,0]);

A1 = double(A1)
B1 = double(B1)

% Question c.
eigA1 = eig(A1) % Unstable at [0,0,0,0]

%Stability for the other Equlibrium Points. There are a total of 4
%equlibrium points. We already checked for the stablity of the (0,0,0,0)
%Point
%0,0,0,0
%pi,0,0,0
%0,pi,0,0
%pi,pi,0,0


A2 = jacobian(X_dot,x);
A2= subs(A2,[theta1,theta2,theta1_dot,theta2_dot],[pi,0,0,0]);
A2= double(A2);

eigA2 = eig(A2)  % Stable at [pi,0,0,0]

A3 = jacobian(X_dot,x);
A3= subs(A3,[theta1,theta2,theta1_dot,theta2_dot],[0,pi,0,0]);
A3= double(A3);

eigA3 = eig(A3) % Unstable at [0,pi,0,0]

A4 = jacobian(X_dot,x);
A4 = subs(A4,[theta1,theta2,theta1_dot,theta2_dot],[pi,pi,0,0]);
A4 = double(A4);

eigA4 = eig(A4) % Unstable [pi,pi,0,0]

% Question d.
rankC = rank(ctrb(A1,B1)) % Controllable because it is full-rank

% Question e.
