function Pout = ForceVel_Main2(K,uke,time,Pk)
%% Parameter Estimation Using fmincon

%% Settings

% Some fmincon settings
MaxItr = 1000000;
MaxFunEvals = 1000000;

% Solve for the current amplitude
I = Pk.thresh+uke*(Pk.sat-Pk.thresh);

% Decimate
Dec = 10;
K = K(1:Dec:end);
time = time(1:Dec:end);
I = I(1:Dec:end);
Ts = mean(diff(time));
L = numel(time);

%% Optimization Settings

% Optimization initialization
% z0 = [Pk.sat;Pk.thresh;Pk.alpha;Pk.beta;...
%     Pk.d11;Pk.d12;Pk.d13;Pk.d14;Pk.d15;Pk.d16;Pk.phik0;...
%     Pk.c23]; 
z0 = [Pk.sat;1;1;Pk.taua];

% Optimization Settings
options = optimset('Algorithm','interior-point',...
    'Display','iter',...
    'MaxFunEvals',MaxFunEvals,'MaxIter',MaxItr,...
    'MaxSQPIter',100); 

%% Linear Constraints and Bounds

% Linear inequality constraints A*x<=B
A = [];
B = [];
% Linear equality constraints Aeq*x=Beq
Aeq = [];
Beq = [];

% Lower bounds
LB = [20;0.001;0.001;0.01];

% Upper bounds
UB = [60;1000;1000;0.5];

%% Parameters for the Optimization

% Store some variabels into a structure to be used by the optimization
optp.L = L;
optp.K = K;
optp.Ts = Ts;
optp.time = time;
optp.x0 = [K(1);0;0];
optp.I = I;
optp.Pk = Pk;

%% Run the Optimization

res = fmincon(@ForceVelCost2,z0,A,B,Aeq,Beq,LB,UB,[],options,optp);

% Extract parameters and package for output
Pout = Pk;
Pout.sat = res(1);
Pout.a4 = res(2);
Pout.a5 = res(3);
Pout.taua = res(4);

% Solve for the zero input joint angle equilibrium
syms peq
Eq = @(peq) Pout.beta/Pout.alpha*cos(peq) - (Pout.d11*(peq-Pout.phik0)+Pout.d13*exp(Pout.d14*peq)-Pout.d15*exp(Pout.d16*peq));
Pout.phieq = double(solve(Eq(peq)==0,peq));


