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



end

%
% STEP #2: Modify, but do NOT rename, this function. It is called every
% time through the simulation loop.
%

function [actuators, data] = runControlSystem(sensors, references, parameters, data)

actuators.l1 = 0;
actuators.l2 = 1;

end
