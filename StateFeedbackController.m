clear; clc; close all;
syms m1 m2 theta1 theta2 r1 r2 l1 l2 I1 I2 theta1_dot theta1_ddot theta2_dot theta2_ddot u1 u2 g 'real'

m1=1; m2=1; l1=1; l2=1 ;r1=0.45; r2=0.45; g=9.81 ;I2= 0.084; I1= 0.084;

% Finding the equilibrium points 
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

% The symbolic reprensatation of the state and input matrix 
u = [u1;u2];
x = [theta1,theta2,theta1_dot,theta2_dot];
x1_dot = theta1_dot;
x2_dot= theta2_dot;
x3_dot= (I2*u1 - I2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) + I2*g*l1*m1*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*l1*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
x4_dot= -(I2*u1 - I1*u2 - I2*u2 - l1^2*m2*u2 - m1*r1^2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1^3*m2^2*r2*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + I2*g*l1*m1*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*u1*cos(theta2) - 2*l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*l1*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + g*l1^2*m1*m2*r2*cos(theta2)*sin(theta1) + l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);

X_dot = [x1_dot; x2_dot; x3_dot; x4_dot];

A = jacobian(X_dot, x);
B = jacobian(X_dot, u);

A = subs(A,[theta1,theta2,theta1_dot,theta2_dot],[0,0,0,0]);
B = subs(B,[theta1,theta2,theta1_dot,theta2_dot],[0,0,0,0]);

A = double(A)
B = double(B)

% Stability around each equilibrium point 
eigA1 = eig(A) % Unstable at [0,0,0,0]

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

% Controllability 
rankC = rank(ctrb(A,B)) % Controllable because it is full-rank

% State-feedback control
syms k11 k12 k13 k14 k21 k22 k23 k24 lambda
K = [k11,k12,k13,k14; k21,k22,k23,k24];
% Acl = A -(B*K);
% AclPoly = simplify(det(Acl - lambda*eye(4)))

lambda = [-1, -2, -1 - 1i, -1 + 1i];
Kvalues = place(A, B, lambda)

% Simulation of the system 
T = 10;
[t,y] = ode45(@ode_link, [0,T], [pi/6,pi/4,0,0]);
K = [23.5850,5.8875,5.1470,2.6108; 5.8875,4.9875,1.5443,0.9770];
for i = 1:size(y,1)
    u(i)= -(K(1,1)*y(i,1) + K(1,2)*y(i,2)+ K(1,3)*y(i,3) + K(1,4)*y(i,4));
    g(i)= -(K(2,1)*y(i,1) + K(2,2)*y(i,2)+ K(2,3)*y(i,3) + K(2,4)*y(i,4));
end

figure(2)
plot(t,y);

figure(1)
subplot(2,2,1);
plot(t,y(:,1),'b');
xlabel('t', 'FontSize',14)
ylabel('theta1','FontSize',14);

subplot(2,2,2);
plot(t,y(:,2),'r');
xlabel('t', 'FontSize',14)
ylabel('theta2','FontSize',14)

subplot(2,2,3);
plot(t,y(:,3),'b');
xlabel('t', 'FontSize',14)
ylabel('theta1 dot','FontSize',14)

subplot(2,2,4);
plot(t,y(:,4),'r');
xlabel('t', 'FontSize',14)
ylabel('theta2 dot','FontSize',14)

figure(3)
subplot(2,1,1);
plot(t,u,'b');
xlabel('t', 'FontSize',14)
ylabel('u1','FontSize',14);

subplot(2,1,2);
plot(t,g,'b');
xlabel('t', 'FontSize',14)
ylabel('u2','FontSize',14);


function dX = ode_link(t,X)
m1=1;m2=1; l1=1; l2=1 ;r1=0.45; r2=0.45; g=9.81 ;I2= 0.084; I1= 0.084;
dX= zeros(4,1);
X=num2cell(X);
[theta1, theta2, theta1_dot, theta2_dot] = deal(X{:});

if abs(theta1)>2*pi
    theta1= mod(theta1,2*pi);

end

if abs(theta2)>2*pi
    theta2= mod(theta2,2*pi);
end

K = [23.585,5.8875,5.1470,2.6108; 5.8875,4.9875,1.5443,0.9770];
U = -K*[theta1;theta2;theta1_dot;theta2_dot];

u1 = [U(1,:)];
u2 = [U(2,:)];

dX(1) = theta1_dot;
dX(2) = theta2_dot;
dX(3) = (I2*u1 - I2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) + I2*g*l1*m1*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*l1*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);

dX(4) = -(I2*u1 - I1*u2 - I2*u2 - l1^2*m2*u2 - m1*r1^2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1^3*m2^2*r2*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + I2*g*l1*m1*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*u1*cos(theta2) - 2*l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*l1*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + g*l1^2*m1*m2*r2*cos(theta2)*sin(theta1) + l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);

end
