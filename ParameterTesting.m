%% Preamble
clear; close all; clc

% Make sure that you are in the right directory
cd('Z:\Labstation\Krysten\Parameter Estimation')

%% Create Threshold/Saturation Stim Signal

% Dynamometer (the piece the attachments go on) should be rotated 20 
% degrees about the vertical axis.

% Make sure on the laptop that amplitude is set at 10 mA.

np = 23; % # pulses (more pulses = better resolution)
tstim = 2;  % Duration of stimulation step [s]
trest = 2;   % Duration of rest between stimulation [s]

ton = (0:tstim+trest:(tstim+trest)*(np-1)) + 5;
tfw = ton + tstim;
time = 0:0.001:tfw(end)+3; % Total time [s]
tf = time(end);

stims = [50 70 90 110 130 150 170 190 210 230 250 270 290 310 330 350 ...
    370 390 410 430 450 470 490]; % pulsewidths [micros]

stim = zeros(1,numel(time)); % Signal for Simulink
for n = 1:numel(time)
    aa = (ton <= time(n)) & (time(n) <= tfw);
    if any(aa)
        stim(n) = stims(aa);
    end
end

% Now run the model and save torque & stim signals

%% Process to Get Threshold/Saturation

tauraw = S1Threshsat15mA{2}.Values.Data;
fs = 1000; % sampling frequency [Hz]
fc = 6; % cutoff frequency [Hz]
[b,a] = butter(4,fc/(fs/2)); % low-pass filter
tau = filtfilt(b,a,tauraw);

ind = find(time == ton(1),1);
biodexstd = std(tau(1:ind)); % st dev of signal before stimulation
biodexm = mean(tau(1:ind)); % mean of signal before stimulation
tau = tau - biodexm; % adjust baseline
biodexm = 0;

% Take the means of the last second of stimulation for all of the
% stimulation pulses.
tavgs = tfw - 0.5;
TauAvg = zeros(size(tavgs));
for n = 1:numel(tavgs)
   
    ind1 = find(tavgs(n)<=time & time<=tfw(n),1,'first');
    ind2 = find(tavgs(n)<=time & time<=tfw(n),1,'last');
    
    TauAvg(n) = mean(tau(ind1:ind2));
    
end

% Find the index of the torque that is 3 times larger than the standard
% deviation of the load cell signal
threshind = find(TauAvg > (3*biodexstd + biodexm),1,'first');
thresh = stims(threshind);

% Find the first stimulation that produces a torque greater than 99% of
% the maximum torque.
ind1 = find(TauAvg > 0.99*max(TauAvg));
sat = stims(ind1(1));

%% Clear Unneccessary Variables

clear answer minst im maxstim stim time mask_enables message_type mname 
clear format_string Biodexstd Biodexm thresh_time thresh_ind LCmax sat_time
clear sat_ind clear h ind1 ind2 I1 I2 Tau1 Tau2 TauDes
clc

%% Passive Stiffness & Mass + Force-Length Params & Activation Constant

% Test 7 different angles

sat = 350;
tf = 10;
time = 0:0.001:tf;
stim = zeros(length(time),1);
for i = 4001:6000
    stim(i) = sat;
end

% Run the model and save torque and stim signals

%% Damping & Inertial Parameters: Pendulum Test

% Use the dynamometer attachment to position the joint.  Dynamometer input 
% is NOT % required for this trial.

tf = 20;
time = 0:0.001:tf;
stim = zeros(length(time),1);

% Run the model and save q, qdot, and qddot

%% Force-velocity Parameter

% This test also does not require dynamometer input.

T = 8; % period [s]
tf = 3*T;
time = 0:0.001:3*tf;
amp = 0.8*sat; % pulsewidth [micros]

for i = 1:length(time)
    stim = amp.*sin(2*pi*time/T);
end

% Run the model and save q, qdot, and qddot

%% Fatigue & Recovery Time Constants

% Intensity [micros]
stims = 0.8*sat;
% Number of warmup pulses
Nw = 10;
% Number of recovery pulses
Nr = 10;
% Time between warmup pulses [s]
dtw = 10;
% Time between recovery pulses [s]
dtr = 10;
% Duration of warmup and recovery pulses [s]
ptw = 1;
% Duration of fatiguing stimulation [s]
tfat = 120;
% Duty cycle
dc = 0.65;
% Time between fatiguing pulses [s]
dtf = (1-dc)*2;
% Number of fatiguing pulses
Nf = tfat/2;
% Duration of fatiguing pulses [s]
ptf = dc*2;
% Final time
tend = 5+Nw*(dtw+ptw)+dtw+tfat+Nr*(dtr+ptw)+dtr;
time = 0:0.001:tend;

tonw = (0:ptw+dtw:(ptw+dtw)*(Nw-1)) + 5; % warmup FES start times
tfw = tonw + ptw; % warmup FES end times

stim = zeros(1,numel(time)); % Signal for Simulink
for n = 1:numel(time)
    aa = (tonw <= time(n)) & (time(n) <= tfw);
    if any(aa)
        stim(n) = stims;
    end
end

tonf = (tfw(end)+10:ptf+dtf:tfw(end)+10+(ptf+dtf)*(Nf-1)); % FES start times
tff = tonf + ptf; % fatigue FES end times
for n = 1:numel(time)
    aa = (tonf <= time(n)) & (time(n) <= tff);
    if any(aa)
        stim(n) = stims;
    end
end

tonr = (115+tfat):ptw+dtr:(115+tfat+(ptw+dtr)*Nr) + 10; % FES start times
tfr = tonr + ptw; % recovery FES end times
for n = 115+tfat+1:numel(time)  
    aa = (tonr <= time(n)) & (time(n) <= tfr);
    if any(aa)
        stim(n) = stims;
    end
end

tf = time(end);

% Run the model and save torque and stim