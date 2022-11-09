
clear all
close all

addpath('D:\tencent\IP\gpops\lib')
addpath('D:\tencent\IP\gpops\nlp')
  
t0 = 0; 
t0max = 0; 
tf=0.8;
tfMin = 0;
tfMax=1;
%%%Œª÷√
x0=0;
xf=0.1;
xmin=0;
xmax=0.1;

y0=0;
yf=0.05;
ymin=0;
ymax=0.05;

q10=0.6*pi;
q1f=pi/2
q1min=0.3*pi;
q1max=0.6*pi;

q20=-0.6*pi;
q2f=-pi;
q2min=-1.2*pi;
q2max=-0.6*pi;

%%%Œ¢∑÷
dx0=0.1;
dxf=-0.02*0;
dxmin=0;
dxmax=20;

dy0=0;
dyf=0.1;
dymin=-2000;
dymax=2000;

dq10=-0.0;
dq1f=-0.0;
dq1min=-pi;
dq1max =20*pi;
ddq1min=-20*pi;
ddq1max=20*pi;

dq20=0;
dq2f=0;
dq2min=-2*pi;
dq2max=2*pi/3;
ddq2min=-pi; 
ddq2max=pi;

u1min = -10000;
u1max =  10000;
u2min = -10000;


u2max =  10000;
u3min = -10000;
u3max =  10000;

u1Guess=0;
u2Guess=0;
u3Guess=0;

xGuess=0;
yGuess=0;
dq1Guess=0;
dq2Guess=0;


iphase=1;
limits(iphase).time.min = [t0 tfMin];
limits(iphase).time.max = [t0 tfMax];
limits(iphase).state.min(1,:) = [x0 xmin xf];
limits(iphase).state.max(1,:) = [x0 xmax xf];
limits(iphase).state.min(2,:) = [y0 ymin yf];
limits(iphase).state.max(2,:) = [y0 ymax yf];
limits(iphase).state.min(3,:) = [q10 q1min q1f];
limits(iphase).state.max(3,:) = [q10 q1max q1f];
limits(iphase).state.min(4,:) = [q20 q2min q2f];
limits(iphase).state.max(4,:) = [q20 q2max q2f];
limits(iphase).state.min(5,:) = [dx0 dxmin dxf];
limits(iphase).state.max(5,:) = [dx0 dxmax dxf];
limits(iphase).state.min(6,:) = [dy0 dymin dyf];
limits(iphase).state.max(6,:) = [dy0 dymax dyf];
limits(iphase).state.min(7,:) = [dq10 dq1min dq1f];
limits(iphase).state.max(7,:) = [dq10 dq1max dq1f];
limits(iphase).state.min(8,:) = [dq20 dq2min dq2f];
limits(iphase).state.max(8,:) = [dq20 dq2max dq2f];

limits(iphase).control.min(1,1) = u1min;
limits(iphase).control.max(1,1) = u1max;
limits(iphase).control.min(2,1) = u2min;
limits(iphase).control.max(2,1) = u2max;
limits(iphase).control.min(3,1) = u3min;
limits(iphase).control.max(3,1) = u3max;

limits(iphase).parameter.min = [];
limits(iphase).parameter.max = [];

limits(iphase).path.min      = [];
limits(iphase).path.max      = [];
limits(iphase).path3.min      = [0];
limits(iphase).path3.max      = [10000];
limits(iphase).event.min     = [];
limits(iphase).event.max     = [];
limits(iphase).duration.min = 1; 
limits(iphase).duration.max = 5;

guess(iphase).time = [0; 1];
guess(iphase).state(:,1) = [x0; xf];
guess(iphase).state(:,2) = [y0; yf];
guess(iphase).state(:,3) = [q10; q1f];
guess(iphase).state(:,4) = [q20; q2f];
guess(iphase).state(:,5) = [dx0; dxf];
guess(iphase).state(:,6) = [dy0; dyf];
guess(iphase).state(:,7) = [dq10; dq1f];
guess(iphase).state(:,8) = [dq20; dq2f];
guess(iphase).control = [0,0,0;0,0,0]*pi/180;
guess(iphase).parameter = [];


%%%% Optimization %%%%
INPUT.optimize.solver = 'ipopt';   %{'ipopt', 'snopt'}
% INPUT.optimize.tol_mesh = 1e-3;
% INPUT.optimize.tol_opt = 1e-6;
% INPUT.optimize.max_mesh_iter = 5;


setup.name = 'Tarjectory';
setup.limits = limits;
setup.guess = guess;
% setup.functions.continuous = @Dynamic_model;
setup.funcs.cost = 'dynamic_end_point';
setup.funcs.dae  = 'Dynamic_model';
% setup.limits = limits;
% setup.guess = guess;
setup.linkages = [];
setup.derivatives = 'finite-difference';
setup.autoscale = 'on';
% setup.mesh.tolerance = 1e-8;
setup.mesh.iteration = 1e6;
% setup.mesh.nodesPerInterval.min = 4;
% setup.mesh.nodesPerInterval.max = 12;
[output,gpopsHistory] = gpops(setup);



% dq1Guess=0;

% x10 = 4.5;
% x1f = 4.5;
% x1min = 0;
% x1max = L;
% x20 = 0;
% x2f = 0;
% x2min = -10*L;
% x2max =  10*L;
% x30 = 0;
% x3f = 2*pi/3;
% x3min = -pi;
% x3max =  pi;
% x40 = 0;
% x4f = 0;
% x4min = -50;
% x4max =  50;
% x50 = pi/4;
% x5f = pi/4;
% x5min =  0;
% x5max =  pi;
% x60 = 0;
% x6f = 0;
% x6min = -50;
% x6max =  50;

% iphase = 1;
%
%
%
% bounds.phase.initialtime.lower = t0;
% bounds.phase.initialtime.upper = t0;
% bounds.phase.finaltime.lower = tfMin;
% bounds.phase.finaltime.upper = tfMax;
% bounds.phase.initialstate.lower = [x0,y0,dq10,dq20];
% bounds.phase.initialstate.upper = [x0,y0,dq10,dq20];
% bounds.phase.state.lower = [xmin, ymin, q1min, q2min];
% bounds.phase.state.upper = [xmax, ymax, q1max, q2max];
% bounds.phase.finalstate.lower = [xmin, ymin, q1min, q2min];
% bounds.phase.finalstate.upper = [xmax, ymax, q1max, q2max];
% bounds.phase.control.lower = [u1min, u2min, u3min];
% bounds.phase.control.upper = [u1min, u2min, u3min]
% bounds.phase.integral.lower = [dxmax, dymax, dq1max, dq2max];
% bounds.phase.integral.upper = [dxmax, dymax, dq1max, dq2max];
%
% % guess.phase.time = [0; 1];
% % guess.phase.state = [[p0; pMax],[q0; qMax]];
% % guess.phase.control = [0; 0];
% % guess.phase.integral = [7.5];
% guess.phase.state = [xGuess, yGuess, dq1Guess, dq2Guess];
% guess.phase.control = [u1Guess, u2Guess, u3Guess];
% guess.phase.time = [t0; tf];
%
% setup.name = 'Launch-Vehicle-Ascent-Problem';
% setup.functions.continuous = @Dynamic_model;
% setup.functions.endpoint = @dynamic_end_point;
% setup.nlp.solver = 'ipopt';
% setup.nlp.snoptoptions.tolerance = 1e-7;
% setup.bounds = bounds;
% setup.guess = guess;
% % setup.auxdata = auxdata;
% setup.derivatives.supplier = 'sparseFD';
% setup.derivatives.derivativelevel = 'second';
% setup.derivatives.dependencies = 'sparseNaN'; % setup.scales.method = 'automatic-bounds'; setup.mesh.method = 'hp-PattersonRao'; setup.mesh.tolerance = 1e-6; setup.method = 'RPM-Differentiation';
% setup.mesh.method = 'hp-PattersonRao';
% setup.mesh.tolerance = 1e-6;
% setup.method = 'RPM-Differentiation';
%
% [output,gpopsHistory] = gpops(setup);


























%
% iphase = 1;
% Bounds on initial and terminal values of time
% limits(iphase).meshPoints = [-1 1];
% limits(iphase).nodesPerInterval = [20];

































% limits(iphase).time.min = [t0 tfmin];
% limits(iphase).time.max = [t0 tfmax];
% limits(iphase).state.min(1,:) = [x0 xmin xf];
% limits(iphase).state.max(1,:) = [x0 xmax xf];
% limits(iphase).state.min(2,:) = [y0 ymin yf];
% limits(iphase).state.max(2,:) = [y0 ymax yf];
% limits(iphase).state.min(3,:) = [dq10 dq1min dq1f];
% limits(iphase).state.max(3,:) = [dq10 dq1max dq1f];
% limits(iphase).state.min(4,:) = [dq20 dq2min dq2f];
% limits(iphase).state.max(4,:) = [dq20 dq2max dq2f];
%
% limits(iphase).control.min(1,:) = u1min;
% limits(iphase).control.max(1,:) = u1max;
% limits(iphase).control.min(2,:) = u2min;
% limits(iphase).control.max(2,:) = u2max;
% limits(iphase).control.min(3,:) = u3min;
% limits(iphase).control.max(3,:) = u3max;
%
% limits(iphase).parameter.min = [];
% limits(iphase).parameter.max = [];
% limits(iphase).path.min      = [];
% limits(iphase).path.max      = [];
% limits(iphase).event.min     = [];
% limits(iphase).event.max     = [];
%
% guess(iphase).time = [0; 100];
% guess(iphase).state(:,1) = [x0; xf];
% guess(iphase).state(:,2) = [y0; yf];
% guess(iphase).state(:,3) = [q10; q1f];
% guess(iphase).state(:,4) = [q20; q2f];
% guess(iphase).control = [0;0;0]*pi/180;
% guess(iphase).parameter = [];
%
%
% %%%% Optimization %%%%
% INPUT.optimize.solver = 'ipopt';   %{'ipopt', 'snopt'}
% INPUT.optimize.tol_mesh = 1e-3;
% INPUT.optimize.tol_opt = 1e-6;
% INPUT.optimize.max_mesh_iter = 5;
%
%
% setup.name = 'Tarjectory';
% setup.limits = limits;
% setup.guess = guess;
% setup.functions.continuous = @Dynamic_model;
% setup.funcs.cost = 'Dynamic_model';
% setup.funcs.dae  = 'robotArmDae';
% setup.linkages = [];
% setup.derivatives = 'finite-difference';
% setup.autoscale = 'off';
% setup.mesh.tolerance = 1e-8;
% setup.mesh.iteration = 200;
% setup.mesh.nodesPerInterval.min = 4;
% setup.mesh.nodesPerInterval.max = 12;
% [output,gpopsHistory] = gpops(setup);