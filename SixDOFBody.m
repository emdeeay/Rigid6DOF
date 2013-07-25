function statedot = SixDOFBody(state,m,II,XYZ,LMN)
% statedot = SixDOFBody(state,m,II,XYZ,LMN)
%
% Author: Mike Abraham
% Last Update: 23 July 2011
%
% Newton-Euler Equations of Motion for a 6 Degree of Freedom Rigid Body
% undergoing general motion in space
%
% state = [x y z psi theta phi u v w p q r]';
% m = mass;
% II = mass moment of inertia matrix in Body-attached frame
% XYZ = Sum of external forces (body frame)
% LMN = Sum of external moments (body frame)
%
% x,y,z are N frame components of the position of the body CG wrt iNertial 
% psi, theta, phi are the standard Yaw, Pitch, Roll euler angles (radians)
% u,v,w are the B frame components of the velocity of CG wrt iNertial
% p,q,r are the B frame components of angular vel of B wrt N

%% Invert inertia matrix
invII = inv(II);

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
ttheta = stheta/ctheta;
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

%% Kinematic Equations

% [xdot; ydot; zdot] = N frame components of velocity of CG wrt N
% [u; v; w] = B frame components of velocity of CG wrt N
%
%  [B Frame Components] = [TBN]*[N frame Components]
%
% Transpose to create TNB...

xyzdot = [TBN']*[u; v; w];

% StrapDown Equations
pqr2ptp = [0    sphi/ctheta   cphi/ctheta;
           0     cphi            -sphi;
           1   sphi*ttheta   cphi*ttheta];
            
ptpdot = pqr2ptp*[p; q; r];
psidot = ptpdot(1);
thetadot = ptpdot(2);
phidot = ptpdot(3);

%% Dynamic Equations

% Force Equation
pqrcross = [0 -r q; r 0 -p; -q p 0];
uvwdot = XYZ/m - pqrcross*[u; v; w];

% Moment taken about CG
pqrdot = invII*(LMN-pqrcross*II*[p; q; r]);

%% Populate State Derivative
statedot = [xyzdot; psidot; thetadot; phidot; uvwdot; pqrdot];
