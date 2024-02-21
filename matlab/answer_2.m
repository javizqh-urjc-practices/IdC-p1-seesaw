close all
clc

syms theta thetadot thetaddot p pdot pddot ms Jb Js r g u;

% Known variables
%Jb = 2.4;
%Js = 0.040;
%ms = 1.25;
%g = 9.81;
%r = 0.05;

% State space representation x = [theta, p, thetadot, pdot];

% Define the dynamic equations
ecuaciones = [(ms*p + Jb)*thetaddot+2*ms*p*pdot*thetadot+ms*g*p*cos(theta)-u==0,(ms+Js/r^2)*pddot-ms*p*thetadot^2+ms*g*sin(theta)==0];

% Solve thetaddot and pddot
[thetaddot, pddot] = solve (ecuaciones,[thetaddot, pddot]);

% Do the derivative to thetaddot and pddot
pdot = diff(pddot);
thetadot = diff(thetaddot);

% State representation
xdot = [thetadot, pdot, thetaddot, pddot]