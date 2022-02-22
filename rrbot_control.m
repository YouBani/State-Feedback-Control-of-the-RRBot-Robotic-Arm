clear; close; clc;
% ROS Setup
rosinit;

j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');

JointStates = rossubscriber('/rrbot/joint_states');

tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);

tau1.Data = 0;
tau2.Data = 0;

send(j1_effort,tau1);
send(j2_effort,tau2);

client = rossvcclient('/gazebo/set_model_configuration');
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};
req.JointPositions = [deg2rad(30), deg2rad(45)];
resp = call(client,req,'Timeout',3);

tic;
t = 0;
i = 1;

while(t < 10)
    t = toc;

    % reading the joint states
    jointData = receive(JointStates);
    x = [jointData.Position(1);jointData.Position(2);
        jointData.Velocity(1);jointData.Velocity(2)];
    k1=[23.5850,5.8875,5.1470,2.6104];
    k2= [5.8875, 4.9875,1.5543,0.9970];

    % State feedback controller 
    tau1.Data = -k1*x;
    tau2.Data = -k2*x;
    
    send(j1_effort,tau1);
    send(j2_effort,tau2);
    
    % Sample data to be plotted 
    g1(i) = jointData.Position(1);
    g2(i) = jointData.Position(2);
    g3(i) = jointData.Velocity(1);
    g4(i) = jointData.Velocity(2);
    
    u1(i) = tau1.Data;
    u2(i) = tau2.Data;
    
    time(i) = t;
    
    i = i + 1;
end

figure(1)
subplot(2,2,1);
plot(time,g1,'b');
xlabel('t', 'FontSize',14)
ylabel('theta1','FontSize',14);

subplot(2,2,2);
plot(time,g2,'r');
xlabel('t', 'FontSize',14)
ylabel('theta2','FontSize',14)

subplot(2,2,3);
plot(time,g3,'b');
xlabel('t', 'FontSize',14)
ylabel('theta1 dot','FontSize',14)

subplot(2,2,4);
plot(time,g4,'r');
xlabel('t', 'FontSize',14)
ylabel('theta2 dot','FontSize',14)

figure(2)
subplot(2,2,1);
plot(time,u1);
xlabel('t');
ylabel('u1');

figure(2)
subplot(2,2,2);
plot(time,u2);
xlabel('t');
ylabel('u2');

tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);
% disconnect from roscore
rosshutdown;