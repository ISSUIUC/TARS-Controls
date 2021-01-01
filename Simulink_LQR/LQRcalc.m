% First load all the variables from LQRsetupvars.mat
% Some of the vars might change when using this script,
% so make sure to update the .mat file when done.
% To run the simulation, run the Simulink file. 

load LQRsetupvars.mat

%Setup variables. These can be changed, but mostly will remain constant.
w = 0.0508;% Width of the flap in m, 2 inches
vi = 279; %Velocity at burnout, in m/s
CDa = 0.5; %CD of rocket aircframe
CD_f = 0.6; %CD of flaps
CL_f = 1.2; %CL of flaps
D = 0.1524; %Diameter of rocket in m, 6 inches
Srefa = pi*D^2/4; %Reference area of the rocket airframe, horizontal crosssection
rho = 1.225; %Constant density of air
m = 45.06; %Mass of rocket after burnout, in kg
I = 0.1552; %Longitudinal moment of inertia, kg m^2
g = 9.81; %Acceleration due to gravity
hd = 3300; %Desired final altitude, in m
hi = 872.21; %Altitude at burnout, in m

%Values for cost function weights and LQR.
%Cost function: int(k1*(h-hd)^2 + k2*(omega)^2 + k3*v^2 + k4*(l1^2 + l2^2))
%where l1 and l2 are the inputs of flap lengths

k1 = 10;
k2 = 1600;
k3 = 700;
k4 = 20;
Q = [k1 0 0; 0 k2 0; 0 0 k3];
R = k4*diag([1 1]);
N = [0 0;0 0;0 0];
[K,S,e] = lqr(linsys1,Q,R,N)

%Additional features that may be changed for simulation purposes
servodelay = 0.08; %Delay to emulate delay in implementation
lmax = 0.0508; %Maximum amount a flap can stretch to

