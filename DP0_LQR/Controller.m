function func = Controller
% INTERFACE
%
%   sensors
%       .t          (time)
%       .z          (position)
%       .zdot       (velocity)
%
%   references
%       .z          (desired position)
%
%   parameters
%       .tStep      (time step)
%       .g          (acceleration of gravity)
%       .m          (mass)
%       .maxthrust  (maximum allowable net thrust; the minimum is always 0)
%
%   data
%       .whatever   (yours to define - put whatever you want into "data")
%
%   actuators
%       .thrust     (net thrust)

% Do not modify this function.
func.init = @initControlSystem;
func.run = @runControlSystem;
end

%
% STEP #1: Modify, but do NOT rename, this function. It is called once,
% before the simulation loop starts.
%

function [data] = initControlSystem(parameters, data)

% Fetch providied scalar parameters
g = parameters.g;
m = parameters.m;
max_act = parameters.max_actuation;
J = parameters.J;
CD_f = parameters.CD_f;
CL_f = parameters.CL_f;
W = parameters.W;
CD = parameters.CD;
rho = parameters.rho;
D = parameters.D;
L = parameters.L;
Sref = parameters.Sref;

% Construct EOMs

syms z zdot phi omega l1 l2

w = [z ; zdot ; phi ; omega];
p = [l1 ; l2];

f_func = [zdot ; ...
            -((rho*zdot^2)/(2*m))*Sref*CD - (rho*zdot^2*CD_f*W*(l1 + l2))/m - g ; ...
            omega ; ...
            (rho*zdot^2*CL_f*D*W*(l1 - l2)) / (2*J) ];
        
h_func = [z ; zdot; omega];

% TODO - Find equilibrium using fsolve instead
state_eqb = [3400 ; 279 ; 0 ; 0.1];
output_eqb = [0.1 ; 0.1];

A = double(subs(jacobian(f_func, w), [w ; p], [state_eqb ; output_eqb]));
B = double(subs(jacobian(f_func, p), [w ; p], [state_eqb ; output_eqb]));
C = double(subs(jacobian(h_func, w), [w ; p], [state_eqb ; output_eqb]));
D = double(subs(jacobian(h_func, p), [w ; p], [state_eqb ; output_eqb]));

k1 = 10;
k2 = 100;
k3 = 500;
k4 = 10;

Qc = [k1 0 0 0; 0 k2 0 0; 0 0 0.0001 0; 0 0 0 k3];
Rc = k4*diag([1 1]);
N = [0 0;0 0;0 0];
K = lqr(A,B,Qc,Rc);

K = 0.00028*K;

data.K = K;

data.state_eqb = state_eqb;
data.output_eqb = output_eqb;

end

%
% STEP #2: Modify, but do NOT rename, this function. It is called every
% time through the simulation loop.
%

function [actuators, data] = runControlSystem(sensors, references, parameters, data)

x = [sensors.z ; sensors.zdot ; 0 ; sensors.omega];

u = -data.K * x

if u(1) > parameters.max_actuation
    u(1) = parameters.max_actuation;
end

if u(1) < 0
    u(1) = 0;
end

if u(2) > parameters.max_actuation
    u(2) = parameters.max_actuation;
end

if u(2) < 0
    u(2) = 0;
end

actuators.l1 = u(1);
actuators.l2 = u(2);

end
