function [ cineq,ceq ] = nonlincon3( z,optp )

%% Nonlinear Inequality Constraints
cineq = [];

%% Extract the Variables

L = optp.L;
x0 = optp.x0;
u = optp.Tau;
par = optp.par;
Ts = optp.Ts;

Kact = optp.Kang;
Kdotact = optp.Kdot;

%% Extract the Variables

% Extract Parameters
par.d12 = z(1);
par.alpha = z(2);
par.beta = 9.81*par.alpha*z(3)*z(4);

%% Evaluate the Dynamics

% Function of the dynamics with stim and motor torque input.
F = @(x,u) [x(2);...
    par.beta*cos(x(1))+par.alpha*(u-(par.d11*(x(1)-par.phik0)+par.d12*x(2)+par.d13*exp(x(1)*par.d14)-par.d15*exp(x(1)*par.d16)))];

x = [x0 zeros(2,L-1)];
for n = 1:L-1;
    x(:,n+1) = RK4Step(F,x(:,n),u(n),Ts);
end

K = x(1,:);
Kdot = x(2,:);

%% Nonlinear Equality Constraints

% Terminal position and velocity constraint
ceq = [Kact(end)-K(end);...
    Kdotact(end)-Kdot(end)];
    

%% End of the Function
end
