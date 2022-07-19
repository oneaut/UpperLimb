%% Fatigue Parameter Estimator

%% Preamble
clear all; close all

addpath('Estimation Extras')

%% Load the Data

% Load the subject data structure
load('SData')

SIDS = fieldnames(SData);
[s,~] = listdlg('PromptString','Select a Subject ID:',...
    'SelectionMode','single',...
    'ListString',SIDS,'ListSize',[160 300]);
SID = char(SIDS(s));

% Location of the results files
resloc = 'Fatigue Results';

evalc(['aa = isfield(SData.' SID  ',''allresFat'');']);
if aa
    evalc(['NF = numel(SData.' SID '.allresFat);']);
else
    error(['No results can be found for ' SID '!'])
end

if NF == 1
    % Only one result file, so load that one
    evalc(['fname = SData.' SID '.allresFat{1};']);
else
    % Ask user to select a results file to load
    evalc(['fnames = SData.' SID '.allresFat;' ]);
    [s,~] = listdlg('PromptString','Select a Results File:',...
        'SelectionMode','single',...
        'ListString',fnames,'ListSize',[200 300]);
    fname = char(fnames(s));
end

% Load the data
load([resloc '\' fname]) 

%% Necessary Data Corrections

% This variable is used to switch the sign on the load cell measurements if
% it is incorrect. 
LCsgn = 1;

% See corrections_fatigue.m in the Estimation Extras directory for the
% corrections that are being made for each data set
corrections_fatigue

%% Prepare the Load Cell Data

% Settings
Ts = 0.001;
dec = 10;

% Trim limits
trimf1 = 129/Ts+1;
trimf2 = (200-15)/Ts+1;
trimr1 = 245/Ts+1;
trimr2 = 320/Ts+1;

% Filter the data
[bf,af] = butter(10,0.1);
mu = filtfilt(bf,af,LCsgn*Tau);

% Remove the bias
mu = mu-mean(mu(1:5001));
mu2 = mu;

% Trim and decimate the data
Ts = dec*Ts;
mu = mu(trimf1:dec:trimf2);
t = time(trimf1:dec:trimf2)-time(trimf1);
a = u(trimf1:dec:trimf2);
% Smooth and normalize
mu2 = mu2/max(mu);
mu = mu/max(mu);
mu = smooth(mu,20);

mur = mu2(trimr1:trimr2);
timeR = time(trimr1:trimr2)-time(trimr1);

L = numel(mu);

%% Muscle Fatigue
% Determine the parameters Tf and mumin from the data collected during the
% fatiguing procedure.

% Time response of muscle fatigue
F = @(p,xdata) p(1)*(1-exp(-xdata/p(2)))+exp(-xdata/p(2));

p0 = [1 2.2];
lb = [0 0];
ub = [1 1000];

p = lsqcurvefit(F,p0,t,mu,lb,ub);

mm = zeros(size(t));
for n = 1:numel(t)
    mm(n) = F(p,t(n));
end

figure(1)
plot(t,mu,t,mm,'--')
axis tight
xlabel('Time [s]')
ylabel('\phi','FontSize',24)
legend('Measured','Fit')
title(['\mu min = ' num2str(p(1)) ',  T_f = ' num2str(p(2))])
grid on

par.mumin = p(1);
par.Tf = p(2);

%% Muscle Recovery
% Determine the parameter Tr from the data collected during the recovery
% procedure.

% Number of points to use for recovery data
np = 6;

[pks,locs] = findpeaks(mur,'minpeakdistance',10/0.001+1);

figure(2)
plot(timeR,mur,timeR(locs),pks,'*')
xlim([min(timeR) max(timeR)])
xlabel('Time [s]')
ylabel('\mu')
title('Recovery Data')
legend('Measured','Peaks','Location','NW')

mudata = pks(1:np);
tdata = timeR(locs(1:np));

tdata = tdata-tdata(1);

% Solution
F = @(p,xdata) 1+(mudata(1)-1)*exp(-xdata/p);

p20 = 30;

lb = 0;
ub = 1000;

p2 = lsqcurvefit(F,p20,tdata,mudata,lb,ub);

mm2 = zeros(size(tdata));
for n = 1:numel(tdata)
    mm2(n) = F(p2,tdata(n));
end

figure(3)
plot(tdata,mudata,'.',tdata,mm2,'MarkerSize',25)
axis tight
legend('Measured','Fit','Location','SE')
title(['T_r = ' num2str(p2)])
xlabel('Time [s]')
ylabel('\phi','FontSize',24)
grid on

par.Tr = p2;

%% Save the Results

choice = questdlg('Would you like to save the data?', ...
    'Save?', ...
    'Yes','No','No');

if strcmpi(choice,'Yes')
    evalc(['SData.' SID '.params' Leg(1) '.Tf = par.Tf;']);
    evalc(['SData.' SID '.params' Leg(1) '.mumin = par.mumin;']);
    evalc(['SData.' SID '.params' Leg(1) '.Tr = par.Tr;']);
    
    save('SData','SData')
end

