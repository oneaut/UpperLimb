function Pout = ForceVel_Main(K,uke,time,Pin)
%% Parameter Estimation Using fmincon

%% Settings

% Some fmincon settings
MaxItr = 1000000;
MaxFunEvals = 1000000;

% Initial condition
x0 = [K(1);0;0];

% Solve for the current amplitude
I = Pin.thresh+uke*(Pin.sat-Pin.thresh);

% Decimate
Dec = 10;
K = K(1:Dec:end);
time = time(1:Dec:end);
I = I(1:Dec:end);
Ts = mean(diff(time));
L = numel(time);


%% Optimization Settings

% Parameters to use as variables of the optimization
% pvar = {'thresh','sat','c23'};
pvar = {'thresh','c23'};

% Optimization initialization
z0 = zeros(numel(pvar),1);
for n = 1:numel(z0)
    z0(n) = getfield(Pin,pvar{n});
end

% Optimization Settings
options = optimset('Algorithm','interior-point',...
    'Display','iter','diffminchange',0,...
    'MaxFunEvals',MaxFunEvals,'MaxIter',MaxItr,...
    'MaxSQPIter',100); 

%% Linear Constraints and Bounds

% Linear inequality constraints A*x<=B
A = [];
B = [];
% Linear equality constraints Aeq*x=Beq
Aeq = [];
Beq = [];

% % Bound known parameters
% PkB = z0(3:end)*[1+V 1-V];

% Lower bounds
% LB = [0;100;0];
LB = [10;0];

% Upper bounds
% UB = [100;100;5];
UB = [40;5];

%% Parameters for the Optimization

% Store some variabels into a structure to be used by the optimization
names = {'L','K','Ts','time','x0','I','Pin','pvar'};

for n = 1:numel(names)
    evalc(['optp.' names{n} ' = ' names{n} ';']);
end

%% Run the Optimization

res = fmincon(@ForceVelCost,z0,A,B,Aeq,Beq,LB,UB,[],options,optp);

%% Output the Resutls

% Extract parameters and package for output
Pout = Pin;
for n = 1:numel(z0)
    evalc(['Pout.' optp.pvar{n} ' = res(n)']);
end

% Solve for the zero input joint angle equilibrium
syms peq
Eq = @(peq) Pout.beta/Pout.alpha*cos(peq) - (Pout.d11*(peq-Pout.phik0)+Pout.d13*exp(Pout.d14*peq)-Pout.d15*exp(Pout.d16*peq));
Pout.phieq = double(solve(Eq(peq)==0,peq));


