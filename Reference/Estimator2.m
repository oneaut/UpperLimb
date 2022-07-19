%% Preamble

clear all
close all
addpath('Estimation Extras')

%% Load the Data

% Load the subject data structure
load('SData')

SIDS = fieldnames(SData);
[s,~] = listdlg('PromptString','Select a Subject ID:',...
    'SelectionMode','single',...
    'ListString',SIDS);
SID = char(SIDS(s));

% Location of the results files
resloc = 'Results';

eval(['fnames = SData.' SID '.allres;']);
[s,~] = listdlg('PromptString','Select a Results File:',...
    'SelectionMode','single',...
    'ListString',fnames,'ListSize',[200 300]);
fname = char(fnames(s));

% Load the data
load([resloc '\' fname])

% Load Popovic parameters for initialization
pop = load('Estimation Extras\AblePopovicParam.mat');

% Determine normalized stimulation amplitude used for isometric contractions
amp = max(uAct);

%% Necessary Data Corrections

% This variable is used to switch the sign on the load cell measurements if
% it is incorrect. 
LCsgn = 1;

% See corrections.m in the Estimation Extras directory for the corrections
% that are being made for each data set
corrections

%% Perform Some Data Manipulation

% Reorder the data
[qPass,ind1] = sort(qPass);
TauPass = LCsgn*TauPass(ind1);
[qIsoCon,ind1] = sort(qIsoCon);
TauIsoCon = LCsgn*TauIsoCon(ind1);

% Find all values that are within +/-delta of any point
delta = 0.15/2;
    
% Average the passive data
L = numel(qPass);
qAvg = NaN*ones(1,L-2);
% TauAvg = [TauPass(1) NaN*ones(1,L-2) TauPass(end)];
TauAvg = NaN*ones(1,L-2);
for n = 2:L-1
    inds = qPass(n)-delta<=qPass & qPass<=qPass(n)+delta;
    inds(end) = 0;
    qAvg(n) = mean(qPass(inds));
    TauAvg(n) = mean(TauPass(inds));
end
% Remove leftover NaNs
TauAvg(isnan(TauAvg)) = [];
qAvg(isnan(qAvg)) = [];
% Put the hyper extension and hyper flexion points back in.
TauAvg = [TauPass(1) TauAvg TauPass(L)];
qAvg = [qPass(1) qAvg qPass(L)];
    
[qP,uind,~] = unique(qAvg);
TauP = TauAvg(uind);

% Average the isometric contraction data
L = numel(qIsoCon);
qAvg = NaN*ones(1,L);
TauAvg = NaN*ones(1,L-2);
for n = 1:L
    inds = qIsoCon(n)-delta<=qIsoCon & qIsoCon<=qIsoCon(n)+delta;
    qAvg(n) = mean(qIsoCon(inds));
    TauAvg(n) = mean(TauIsoCon(inds));
end
    
[qI,uind,~] = unique(qAvg);
TauI = TauAvg(uind);

%% Initialize Parameters Structure
% Some parameters need to be initialized.

% Subject mass and height
eval(['H = SData.' SID '.height;']);
eval(['M = SData.' SID '.weight;']);
if ~isnumeric(H)
    H = str2double(H);
end
if ~isnumeric(M)
    M = str2double(M);
end

% Gravitational acceleration
g = 9.81;

par.sat = sat;
par.thresh = thresh;
par.taua = 0.01;

% Initialize muscle model with Popovic1999 parameters
par.d11 = pop.d11;
par.d12 = pop.d12;  
par.d13 = pop.d13;  
par.d14 = pop.d14;  
par.d15 = pop.d15;  
par.d16 = pop.d16;
par.phik0 = pop.phik0;
par.c23 = pop.c23;
par.c22 = pop.c22;
par.c21 = pop.c21;
par.c20 = pop.c20;
par.a1 = 248.7;
par.a2 = 1.382;
par.a3 = 0.820;
par.a4 = 0.5868;
par.a5 = 1.1054;

% Initialize using anthropometric data
par.lc = 0.1727*H;
par.m = 0.061*M;
par.J = 8.57E-4*M*H^2;
par.alpha = 1/par.J;
pop.alpha = par.alpha;
par.beta = par.alpha*par.m*g*par.lc;
pop.beta = par.beta;

% Structure of initialized parameters
pinit = par;

%% Passive Data Using Anthropometry

% If the hyperextension passive data point was negative make it equal to
% zero.
% qP(qP<0) = 0;

% Fit the passive data
% Initial guess
p0 = [par.d11 par.phik0 par.d13 par.d14 par.d15 par.d16];
p0 = zeros(size(p0));

% Bounds on the parameters
lb = [0 0 0 0 0 -inf];
ub = [inf inf inf inf inf 0];

% Function of the stiffness and mass dynamics
F = @(p,x) par.beta/par.alpha*cos(x) - (p(1)*(x-p(2))+p(3)*exp(p(4)*x)-p(5)*exp(p(6)*x));

% Use the lsqcurvefit function to determine the gravitational and 
options=optimset('TolFun',1e-12,'TolX',1e-12,'MaxFunEvals',40000,'MaxIter',50000);
[p,~,~,~,~] = lsqcurvefit(F,p0,qP,TauP,lb,ub,options);

% Evaluate and plot the results
phiP = linspace(min(qP)-0.05,max(qP)+0.05,50);
TresP = F(p,phiP);

% Update the parameters structure
par.d11 = p(1); 
par.phik0 = p(2);
par.d13 = p(3);
par.d14 = p(4);
par.d15 = p(5);
par.d16 = p(6);

% Plot the results
figure(1)
plot(phiP*180/pi,TresP,qP*180/pi,TauP,'*',qPass*180/pi,TauPass,'o','MarkerSize',10)
legend('Fit','Averaged','Measured','Location','Best')
axis tight

%% Check Saturation and Threshold
% Plot the stim ramp data to check the saturation and threshold levels that
% were already computed

% Preliminary calculations
ind = 3001;
rampTau = LCsgn*rampTau;
LCstd = std(rampTau(1:ind));
rTau = rampTau-mean(rampTau(1:ind));
%%
figure(11)
[ax,p1,p2] = plotyy(ramptime,rTau,...
    ramptime,reshape(rampstim/1000,1,numel(ramptime)));
ylabel(ax(1),'Isometric Torque [Nm]','FontSize',24)
ylabel(ax(2),'Current [mA]','FontSize',24)
xlabel('Time [s]','FontSize',24)
ylim(ax(1),[0 1.05*max(rTau)])
ylim(ax(2),[0 1.05*max(rampstim)/1000])
xlim(ax(1),[0 max(ramptime)])
xlim(ax(2),[0 max(ramptime)])
set(p2,'LineStyle','--')
grid on

%% Force-Length Parameters

F = @(p,xdata) amp*(p(1)*exp(-(xdata-p(2)).^2/(2*p(3))));

% Make an initial guess at the parameters
% Paramter array: [a1 a2 a3]
p0 = [max(-TauI) pinit.a2 pinit.a3];

% Bounds on the parameters
lb = [0.001 0.001 0.001];
ub = [inf inf inf];

% Options for lsqcurvefit
options = optimset('MaxFunEvals',1000);

% Use the lsqcurvefit function to determine the gravitational and 
[p,~,~,~,~] = lsqcurvefit(F,p0,qI,-TauI,lb,ub,options);

% Store the results in the structure
par.a1 = p(1);
par.a2 = p(2);
par.a3 = p(3);

% Evaluate and plot the results
phiI = linspace(0,max(qI)*1.1*180/pi,100)*(pi/180);
TresI = F(p,phiI);
%%
% Plot the curve fit
figure(2)
plot(phiI*(180/pi),TresI,qI*(180/pi),-TauI,'*',qIsoCon*(180/pi),-TauIsoCon,'o','MarkerSize',10)
xlabel('Joint Angle [deg.]','FontSize',24)
ylabel('Torque [Nm]','FontSize',24)
legend('Fit','Averaged','Measured','Location','Best')
axis tight

%% Muscle Activation Time Constant

% Select a trial number
tnum = 5;

% Select times to trim the data
% New trim settings
trim1 = 2001;
trim2 = 5001;

switch SID
    case 'S1'
        switch fname
            case 'Results_Mar_30_2015_14_29'
                % Good trial number for Nick's 3/30 data
                tnum = 15;
                trim1 = 5001;
                trim2 = 8001;
            case 'Results_May_11_2015_15_49'
                % Good trial number for Nick's 3/30 data
                tnum = 15;
        end
        
    case 'S2'
        switch fname
            case 'Results_Apr_08_2015_17_29'
                % Good trial number for Naji's 4/8 Data
                tnum = 14;
                trim1 = 5001;
                trim2 = 8001;
            case 'Results_Apr_14_2015_16_35'
                % Good trial number for Naji's 4/16 Data
                tnum = 18;
                trim1 = 5001;
                trim2 = 8001;
            otherwise
                tnum = 1;
        end
    case 'S4'
        % Good trial number for Marcus' 5/5 data
        tnum = 11;
    case 'S6'
        % Good trial number for Ricardo's 5/4 data
        tnum = 10;
    otherwise
        tnum = 1;
end

% Remove the LC bias and trim the data
eval(['T = LCsgn*TauIsoCon' num2str(tnum) '(trim1:trim2)-mean(TauIsoCon' num2str(tnum) '(1:trim1));']);
% Smooth the data
T = smooth(T,20);
% Lowpass filter the load cell data
[b,a] = butter(10,0.15);
T = -filtfilt(b,a,T);
eval(['uke = uIsoCon' num2str(tnum) '(trim1:trim2);']);
% Remove the bias
T = T-mean(T(1:1001));
ake = (T/max(T(1001:end)))*amp;
tact = timeAct(trim1:trim2)-timeAct(trim1);

% Create a range of time constant values to check for a minimum
tauavec = 0:0.01:1;
err = zeros(size(tauavec));
for n = 1:numel(tauavec)  
    H = tf(1,[tauavec(n) 1]);
    [yact,~,~] = lsim(H,uke,tact);
    err(n) = rms(ake-yact);
end

[~,I] = min(err);
%%
% Plot the error vs. the activation time constant
figure(3)
plot(tauavec,err,tauavec(I),err(I),'v','MarkerSize',10)
xlabel('Activation Time Constant [s]','FontSize',24)
ylabel('RMS Error','FontSize',24)
title(['\tau_a = ' num2str(tauavec(I))])

figure(4)
H = tf(1,[tauavec(I) 1]);
[yact,~,~] = lsim(H,uke,tact);
plot(tact,ake,tact,yact,'--')
title(['\tau_a = ' num2str(tauavec(I))])
axis tight
xlabel('Time [s]','FontSize',24)
ylabel('Normalized Activation','FontSize',24)
legend('Measured','Fit','Location','SE')
grid on

% Save the computed taua parmaeter
par.taua = tauavec(I);

%% Pendulum Test Using fmincon

% Edit the data
dec = 10;
Ts = 0.001*dec;
K = qPend(1:dec:end);
Tau = TauPend(1:dec:end);
% Lowpass filter the load cell and joint angle data
[b,a] = butter(10,0.15);
K = filtfilt(b,a,K);
Tau = filtfilt(b,a,Tau);
% Smooth the data
K = smooth(K,20);
Tau = smooth(Tau,40);
Kdot = smooth(diff([K(1);K])/Ts);
% Trim the data
ind1 = find(abs(Kdot) > 0.2,1,'first');
ind2 = find(abs(Kdot) > 0.01,1,'last');
% ind2 = length(K);
K = K(ind1:ind2);
Kdot = Kdot(ind1:ind2);  
Tau = Tau(ind1:ind2);
tPend = 0:Ts:(numel(K)-1)*Ts;

for n = 1:1
    [d12,alpha,m,lc] = KneePendEstimation_Main3(K,Kdot,Tau,tPend,par);
end

% Evaulate the dynamics with the new parameters
% Function of the dynamics with only motor torque input.
F = @(x,u) [x(2);...
    9.81*alpha*m*lc*cos(x(1))+alpha*(u-(par.d11*(x(1)-par.phik0)+d12*x(2)+par.d13*exp(x(1)*par.d14)-par.d15*exp(x(1)*par.d16)))];
x0 = [K(1);Kdot(1)];
L = numel(K);
x = [x0 zeros(2,L-1)];
for n = 1:L-1;
    x(:,n+1) = RK4Step(F,x(:,n),Tau(n),Ts);
end

Kfit = x(1,:);
%%
figure(5)
% plot(tPend,[K Kfit']*180/pi)
hold on
plot(tPend,K*180/pi,'Linewidth',2)
plot(tPend,Kfit'*180/pi,'g--','Linewidth',2)
hold off
xlabel('Time [s]','FontSize',24); ylabel('Joint Angle [deg.]','FontSize',24)
legend('Measured','Fit','Location','SE')
axis tight

% Update the parameters
par.alpha = alpha;
par.d12 = d12;
par.m = m;
par.lc = lc;
par.beta = 9.81*alpha*m*lc;

%% Adjust the Mass and Inertial Parameters

% Adjust the mass and moment of inertia parmeters
% Mass of the leg extension machine [kg]
ml = 5.24;
% Length to center of mass of the leg extension machine [m]
lcl = 0.266;
% Moment of inertia of the leg extension machine [kgm^2]
Jl = 0.558;
% Adjusted mass
mnew = ml+par.m;
% Adjusted length to center of mass
lcnew = (par.lc*par.m+lcl*ml)/mnew;
% Adjusted moment of inertia
Jnew = par.J+Jl+par.m*(par.lc^2-lcnew^2)+ml*(lcl^2-lcnew^2);

% Compute the new alpha and beta parameters
par.alpha = 1/Jnew;
par.beta = mnew*g*lcnew*par.alpha;

% Solve for the equilibrium
syms peq
Eq = @(peq) par.beta/par.alpha*cos(peq) - (par.d11*(peq-par.phik0)+par.d13*exp(par.d14*peq)-par.d15*exp(par.d16*peq));
par.phieq = double(solve(Eq(peq)==0,peq));

%% Force-Velocity Parameter Estimation

% Trim the data
trim = 30001;
uS = uSine(1:trim)';
Tau = -TauSine(1:trim)';
qS = qSine(1:trim);
qS = qS-qS(1)+par.phieq;
L = numel(uS);
Ts = 0.001;
x0 = [qS(1);0;0];
x = [x0 zeros(3,L-1)];
tS = 0:Ts:(L-1)*Ts;

% Use data to adjust the threshold first
% Solve for the current
I = pinit.thresh+uS*(pinit.sat-pinit.thresh);
ind = find(qS*180/pi<(qS(1)*180/pi-0.25),1);
par.thresh = I(ind);
% Recompute u
uS = (I-par.thresh)/(par.sat-par.thresh);

% Use an optimization to compute the best sat, thresh, and c23 parmeters
Pout = ForceVel_Main2(qS,uS,tS,par);

% Resolve for the normalized stimulation using the new saturation and
% threshold.
I = par.thresh+uS*(par.sat-par.thresh);
u = (I-Pout.thresh)/(Pout.sat-Pout.thresh);
u(u<0) = 0; u(u>1) = 1;

% Solve for dynamics using the best c23 parameter
% Redetermine the dynamics function
F = @(x,u) [x(2);...
    Pout.beta*cos(x(1))-Pout.alpha*((Pout.a1*exp(-(x(1)-Pout.a2).^2/(2*Pout.a3)))*(Pout.a4/(1+(Pout.a4-1)*exp(-Pout.a5*x(2))))*x(3)+...
    (Pout.d11*(x(1)-Pout.phik0)+Pout.d12*x(2)+Pout.d13*exp(x(1)*Pout.d14)-Pout.d15*exp(x(1)*Pout.d16)));...
    (u-x(3))/Pout.taua];

for m = 1:L-1;
    x(:,m+1) = RK4Step(F,x(:,m),u(:,m),Ts);
end
%%
% Plot the results
figure(10)
plot(tS,qS'*180/pi,tS,x(1,:)*180/pi,'--')
set(gca,'YDir','reverse');
axis tight
xlabel('Time [s]','FontSize',24)
ylabel('Knee Joint Angle [deg.]','FontSize',24)
legend('Measured','Model','Location','SE')
erms = rms((x(1,:) - qS')*180/pi);
grid on
% title(['a4: ' num2str(Pout.a4) ', a5: ' num2str(Pout.a5) ', RMS: ' num2str(erms) ' degrees'])
%%
% Update the par array
par.thresh = Pout.thresh;
par.sat = Pout.sat;
par.a4 = Pout.a4;
par.a5 = Pout.a5;
par.taua = Pout.taua;

%% Store the Estimated Parameters

choice = questdlg('Would you like to store and save the results into SData?', ...
	'Save Results', ...
	'Yes','No','No');

if strcmp(choice,'Yes')
    % Names of parameters to add
    pNames = {'a1','a2','a3','a4','a5'};
    for n = 1:numel(pNames)
       eval(['SData.' SID '.params' Leg(1) '.' pNames{n}...
           ' = par.' pNames{n} ';']);
    end
    % Save the SData structure
    save('SData','SData')
end
