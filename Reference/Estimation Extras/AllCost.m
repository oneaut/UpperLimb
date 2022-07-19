function cost = AllCost( z,optp )

%% Extract the Variables

time = optp.time;
Ts = optp.Ts;
x0s = optp.x0s;
x0p = optp.x0p;
I = optp.I;
Pk = optp.Pk;

% Push/pull data
qPP = optp.qPP;
TauPPact = optp.TauPP;

% Actual knee angles
Kpact = optp.Kp;
Ksact = optp.Ks;

%% Extract the Variables

% Extract parameters that are being adjusted
P.sat = z(1);
P.thresh = z(2);
P.alpha = z(3);
P.beta = z(4);
P.d11 = z(5);
P.d12 = z(6);
P.d13 = z(7);
P.d14 = z(8);
P.d15 = z(9);
P.d16 = z(10);
P.phik0 = z(11);
P.c23 = z(12);

% Parameters that are not being adjusted
P.c22 =Pk.c22;
P.c21 = Pk.c21;
P.c20 = Pk.c20;
P.taua = Pk.taua;

%% Push/Pull Data

TauPP = P.beta/P.alpha*cos(qPP) - (P.d11*(qPP-P.phik0)+P.d13*exp(P.d14*qPP)-P.d15*exp(P.d16*qPP));

%% Pendulum Dynamics

options = odeset();
[~,xp] = ode45(@penddyn,time,x0p,options,P);

Kp = xp(:,1);

%% Sinusoidal Stimulation Dynamics

% Function of the dynamics
F = @(x,u) [x(2);...
    P.beta*cos(x(1))-P.alpha*((P.c22*x(1)^2+P.c21*x(1)+P.c20)*(1+P.c23*x(2))*x(3)+...
    (P.d11*(x(1)-P.phik0)+P.d12*x(2)+P.d13*exp(x(1)*P.d14)-P.d15*exp(x(1)*P.d16)));...
    (u-x(3))/P.taua];

% Comput u using the update saturation and threshold
u = (I-P.thresh)/(P.sat-P.thresh);

Ls = numel(u)-1;
xs = [x0s zeros(3,Ls)];
for n = 1:Ls;
   xs(:,n+1) = RK4Step(F,xs(:,n),u(n),Ts);
end
    
% Extrac the knee joint angle
Ks = xs(1,:);

%% Evaluate the Cost Function

% Quadratic cost function
cost = sum((Kpact-Kp).^2)/numel(Kp) + sum((Ksact-Ks').^2)/numel(Ks) + 10*sum(TauPPact-TauPP).^2/numel(TauPP);

end

