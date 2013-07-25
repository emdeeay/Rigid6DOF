%% Clear the workspace variables, close all figures, clear the command window
clear
close all
clc

%% Initial conditions
xyz0 = [0; 0; 0];
ptp0 = pi/180*[0; 45; 0];
uvw0 = [50; 0; 0];
pqr0 = [.4; .4; .05;];

state0 = [xyz0; ptp0; uvw0; pqr0];

% Timestep
dt = .1;

% Total simulation time
tfinal = 120;

% Construct an array of time values
t = 0:dt:tfinal;

% Construct an array that will store the state history, initialize to zeros
state = zeros(numel(state0),numel(t));

% Place the initial conditions in the beginning of the state array:
state(:,1) = state0;

% The integration loop
for idx = 2:numel(t)
    state(:,idx) = RK4Integrator(@A_Rock,t(idx-1),state(:,idx-1),dt);
end

% unwrap state vector
xyz = state(1:3,:);
ptp = state(4:6,:);
uvw = state(7:9,:);
pqr = state(10:12,:);

% Plots
figure(1)
plot3(xyz(1,:),xyz(2,:),xyz(3,:))
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on

figure(2)
plot(t,ptp)
legend('yaw','pitch','roll')
xlabel('time')
grid on

figure(3)
plot(t,uvw)
legend('u','v','w')
grid on
xlabel('time')

figure(4)
plot(t,pqr)
legend('p','q','r')
grid on
xlabel('time')