%% Muscle Fatigue Parameter Estimation Tests

%% Preamble
close all; clc

% Make sure that you are in the right directory
% cd('D:\Dropbox\Nick\Muscle Parameter Estimation')
addpath('Extras') 

load SData

% Stimulation settings
PW = 400;           % Pulse width of the stimulation [microseconds]
stimfreq = 35;      % Frequency of the stimulation [Hz]
LCconv = 23.2;     % Conversion factor from LC voltage to torque [Nm/V]

%% Subject Information

SID = selectSubject;

Leg = '';
while isempty(Leg)
    Leg = questdlg('Which leg will you run the procedure on?', ...
        'Choose Leg', ...
        'Left','Right','Right');
end

%% Muscle Fatigue Test

Ready(8)

evalc(['sat = SData.' SID '.params' Leg(1) '.sat;']);
evalc(['thresh = SData.' SID '.params' Leg(1) '.thresh;']);

% Number of warmup pulses
Nw = 10;
% Duration of fatiguing stimulation [s]
tfat = 120;
% Number of recovery pulses
Nr = 10;
% Time between warmup pulses [s]
dtw = 10;
% Time between recovery pulses [s]
dtr = 10;
% Duration of warmup and recovry pulses
ptw = 1;

tend = 5+Nw*(dtw+ptw)+dtw+tfat+Nr*(dtr+ptw)+dtr;

%%%%%%%%%%%%%%%%%
mname = 'Simulink Models/FatigueTesting';
% Open the Simulink model
open_system(mname)
% Build the Simulink model
qc_build_model
% Run the Simulink model
qc_connect_model
pause(3)
qc_start_model
% Pause and wait for the model to finish
pause(tend+3)
% Close the Simulink model
save_system
close_system
%%%%%%%%%%%%%%%%%

%% Save Collected Data
% Organize the data that has been collected in the above procedures and
% save it for processing.

stamp = datestr(clock,'yyyy_mm_dd_HH_MM');
fname = ['Results_' stamp];

load SData

evalc(['aa = isfield(SData.' SID  ',''allresFat'');']);
if aa
    evalc(['NF = numel(SData.' SID '.allresFat);']);
else
    NF = 0;
end
evalc(['SData.' SID '.allresFat{NF+1} = fname;']);

% Save the SData structure
save('SData','SData')

clear SData aa 

% Save the resulting data
save(['Fatigue Results\' fname])


