function statedot = A_Rock(t,state)

% Equations of motion for a 6 DOF rock in a gravitational field.
% Simple damping models for linear and rotational velocities are included.

%% Rigid Body Properties
m = 1;
Ixx = 1;
Iyy = 5;
Izz = 10;
Ixy = 0;
Ixz = 0;
Iyz = 0;

%% Construct Inertia Matrix
II = [Ixx Ixy Ixz;
      Ixy Iyy  Iyz;
      Ixz Iyz Izz];
  
%% Unwrap state vector

% Positions wrt iNertial space (use your favorite unit)
xposn = state(1);
yposn = state(2);
zposn = state(3);

% Euler Angles (radians)
psi = state(4);
theta = state(5);
phi = state(6);

% Body Velocities (same length unit as posn, same time unit as angular vel)
u = state(7);
v = state(8);
w = state(9);

% Body Angular Rates (radians/time)
p = state(10);
q = state(11);
r = state(12);

%% Transcendental Functions of Euler Angles
cpsi = cos(psi);
spsi = sin(psi);
ctheta = cos(theta);
stheta = sin(theta);
cphi = cos(phi);
sphi = sin(phi);

%% Populate Rotation Matrices

L3psi = [ +cpsi +spsi   0; ...
          -spsi +cpsi   0; ...
            0     0     1];

L2theta = [ +ctheta   0   -stheta; ...
                0     1      0   ; ...
            +stheta   0   +ctheta];

L1phi = [  1    0     0   ; ...
           0  +cphi  +sphi; ...
           0  -sphi  +cphi];

% Transformation from iNertial frame (N) to body frame (B)
TBN = L1phi*L2theta*L3psi;       
 
%% External Forces and Moments

% The gravity vector has magnitude g and points in the Down direction.
% Therefore, it is easy to write in the N frame:
Gravity_N = [0; 0; 9.81];

% Transform into the Body frame:
Gravity_B = TBN*Gravity_N;

% Weight is a force. Gravity is an acceleration
Weight_B = m*Gravity_B;

% Drag is proportional to (and opposite) the velocity vector.
Drag_B = -0.1*[u; v; w];

% Sum of external forces in the body frame
XYZ = Weight_B + Drag_B;

% Sum of external moments

RotDrag_B = -0.02*[p; q; r];

LMN = RotDrag_B;

%% 6DOF Equations
statedot = SixDOFBody(state,m,II,XYZ,LMN);

