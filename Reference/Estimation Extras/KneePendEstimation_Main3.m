function [d12,alpha,m,lc] = KneePendEstimation_Main3(K,Kdot,Tau,time,par)
%% Parameter Estimation Using fmincon

%% Settings

% Some fmincon settings
MaxItr = 1000;

% Decimate
Dec = 1;
K = K(1:Dec:end);
time = time(1:Dec:end);
Tau = Tau(1:Dec:end);
Ts = time(2)-time(1);

%% Optimization Settings

% Optimization initialization
z0 = [0;par.alpha;par.m;par.lc];

% Optimization Settings
options = optimset('Algorithm','interior-point',...
    'Display','iter','MaxFunEvals',10000000,'MaxIter',MaxItr,...
    'MaxSQPIter',100); 

%% Linear Constraints and Bounds

% Linear inequality constraints A*x<=B
A = [];
B = [];
% Linear equality constraints Aeq*x=Beq
Aeq = [];
Beq = [];

% Select how much you will allow mass and length to center of mass to vary
% from their original values.
V = 0.10;
mbnd = par.m*[1+V 1-V];
lcbnd = par.lc*[1+V 1-V];

LB = [0 0 min(mbnd) min(lcbnd)];
UB = [inf inf max(mbnd) max(lcbnd)];

%% Parameters for the Optimization

% Store some variabels into a structure to be used by the optimization
optp.L = numel(K);
optp.Kang = K;
optp.Kdot = Kdot;
optp.time = time;
optp.x0 = [K(1);Kdot(1)];
optp.Tau = Tau;
optp.par = par;
optp.Ts = Ts;

%% Run the Optimization

% res = fmincon(@objcost3,z0,A,B,Aeq,Beq,LB,UB,@nonlincon3,options,optp);
res = fmincon(@objcost3,z0,A,B,Aeq,Beq,LB,UB,[],options,optp);

% Extract the results
d12 = res(1);
alpha = res(2);
m = res(3);
lc = res(4);




