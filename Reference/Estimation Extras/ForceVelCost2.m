function cost = ForceVelCost2( z,optp )

%% Extract the Variables

L = optp.L;
Ts = optp.Ts;
x0 = optp.x0;
I = optp.I;
Pk = optp.Pk;

% Actual knee angles
Kact = optp.K;

%% Extract the Parameters

% Known parameters
P = Pk;
% Update the parameters that are the variables in the optimization
P.sat = z(1);
P.a4 = z(2);
P.a5 = z(3);
P.taua = z(4);

%% Function of the Dynamics

F = @(x,u) [x(2);...
    P.beta*cos(x(1))-P.alpha*((P.a1*exp(-(x(1)-P.a2).^2/(2*P.a3)))*(P.a4/(1+(P.a4-1)*exp(-P.a5*x(2))))*x(3)+...
    (P.d11*(x(1)-P.phik0)+P.d12*x(2)+P.d13*exp(x(1)*P.d14)-P.d15*exp(x(1)*P.d16)));...
    (u-x(3))/P.taua];

%% Evaluate the Dynamics
% The nonlinear equality constraints are the dynamics, which are solved
% using a 4th order Runge-Kutta method.

% Compute u using the update saturation and threshold
u = (I-P.thresh)/(P.sat-P.thresh);
u(u<0) = 0; u(u>1) = 1;

x = [x0 zeros(3,L-1)];
for n = 1:L-1;
   x(:,n+1) = RK4Step(F,x(:,n),u(n),Ts);
end

% Extract the knee joint angle
K = x(1,:);

% Solve for the zero input joint angle equilibrium
syms peq
Eq = @(peq) P.beta/P.alpha*cos(peq) - (P.d11*(peq-P.phik0)+P.d13*exp(P.d14*peq)-P.d15*exp(P.d16*peq));
phieq = double(solve(Eq(peq)==0,peq));

%% Evaluate the Cost Function

% Quadratic cost function
cost = sum((Kact-K').^2)/L + (phieq-Kact(1)).^2;

end