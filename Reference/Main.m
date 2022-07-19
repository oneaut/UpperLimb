%% Parameter Estimation Tests

%% Preamble
clear all; close all; clc

% Make sure that you are in the right directory
cd('Z:\Labstation\Qiang\Muscle Parameter Estimation')
addpath('Extras') 

% Stimulation settings
CR = 25;           % Current of the stimulation [mA]
stimfreq = 35;      % Frequency of the stimulation [Hz]
LCconv = 26.89;     % Conversion factor from LC voltage to torque [Nm/V]

%% Subject Information

SID = selectSubject_New;

Leg = '';
while isempty(Leg)
    Leg = questdlg('Which leg will you run the procedure on?', ...
        'Choose Leg', ...
        'Left','Right','Right');
end
                   
%% Procedure 1) Setup the Encoder

Ready(1)

board_handle = hil_open('qpid_e');

pause(3)

% Instruct the experimenter to move the leg to fully extended
figure(1)
imshow('Extras\Im1.png','Border','Tight')

% Pause then set the encoder value to zero
pause(3)
hil_set_encoder_counts(board_handle,0,0)

% Instruct the experimenter to release the leg
close(1)
figure(2)
imshow('Extras\Im5.png','Border','Tight');
pause(3)
close(2)

% Close the board
hil_close(board_handle);

% Clear any unneccessary variables
clear board_handle

%% Procedure 2) Pendulum Test

Ready(2)

%%%%%%%%%%%%%%%%%
mname = 'Simulink Models/PendulumTest';
% Open the Simulink model
open_system(mname)
% Build the Simulink model
qc_build_model
% Run the Simulink model
qc_connect_model
pause(3)
qc_start_model
% Pause and wait for the model to finish
pause(13)
% Close the Simulink model
save_system
close_system
%%%%%%%%%%%%%%%%%

clear mname

%% Procedure 3) Push/Pull Test

Ready(3)

%%%%%%%%%%%%%%%%%
mname = 'Simulink Models/PushPull';
% Open the Simulink model
open_system(mname)
% Build the Simulink model
qc_build_model
% Run the Simulink model
qc_connect_model
pause(3)
qc_start_model
% Pause and wait for the model to finish
pause(13)
% Close the Simulink model
save_system
close_system
%%%%%%%%%%%%%%%%%

% Clear unnecessary variables
clear mname

%% Procedure 4) Warmup
% This procedure issues pulses to the muscle to warm it up to the
% stimulation.

Ready(4)

%%%%%%%%%%%%%%%%%
mname = 'Simulink Models/Warmup';
% Open the Simulink model
open_system(mname)
% Build the Simulink model
qc_build_model
% Run the Simulink model
qc_connect_model
pause(3)
qc_start_model
% Pause and wait for the model to finish
pause(8.5+max(np)*(5.5))
% Close the Simulink model
save_system
close_system
%%%%%%%%%%%%%%%%%

clear mname format_string mask_enables np amp message_type

%% Procedure 5) Stimulation Ramp

Ready(5)

% Prompt the user for a range of stimulation values to search in.
answer = inputdlg({'Lower Range [ma]:','Upper Range [ma]:'},...
    'Stimulation Range',1,{'20','70'});
minstim = str2double(answer{1});
maxstim = str2double(answer{2});

% Create the stimulation train
% Number of pulses that you wish to use. More pulses gives a better
% resolution for the determination of the threshold and saturations.
np = 30;
tstim = 1;  % Duration of stimulation step [s]
trest = 4;   % Duration of rest between stimulation [s]

ton = (0:tstim+trest:(tstim+trest)*(np-1))+5;
toff = ton+tstim;
m = (maxstim-minstim)/(ton(end)-ton(1));
time = 0:0.001:toff(end)+3;

stims = (ton-5)*m+minstim;

stim = zeros(1,numel(time));
for n = 1:numel(time)
    aa = (ton <= time(n)) & (time(n) <= toff);
    if any(aa)
        stim(n) = stims(aa);
    end
end

%%%%%%%%%%%%%%%%%
mname = 'Simulink Models/StimRamp_V2';
% Open the Simulink model
open_system(mname)
% Build the Simulink model
qc_build_model
% Run the Simulink model
qc_connect_model
pause(3)
qc_start_model
% Pause and wait for the model to finish
pause(time(end)+3)
% Close the Simulink model
save_system
close_system
%%%%%%%%%%%%%%%%%

% Script for processing the data from StimRamp_V2.slx
% Standard deviation of the LC signal noise.
ind = find(time == ton(1),1);
LCstd = std(rampTau(1:ind));
% Mean of the LC signal before stimulation is applied.
LCm = mean(rampTau(1:ind));
% Take the means of the last second of stimulation for all of the
% stimulation pulses.
tavgs = toff-0.5;
TauAvg = zeros(size(tavgs));
for n = 1:numel(tavgs)
   
    ind1 = find(tavgs(n)<=time & time<=toff(n),1,'first');
    ind2 = find(tavgs(n)<=time & time<=toff(n),1,'last');
    
    TauAvg(n) = mean(rampTau(ind1:ind2));
    
end
% Find the index of the torque that is 3 times larger than the standard
% deviation of the load cell signal
threshind = find(TauAvg>(3*LCstd+LCm),1,'first');
thresh = stims(threshind);

% Find the lowest stimulation that produces a torque greater than 95% of
% the maximum torque.
ind1 = find(TauAvg < 100,1,'last');
ind2 = find(TauAvg > 100,1,'first');
I1 = stims(ind1);
I2 = stims(ind2);
Tau1 = TauAvg(ind1);
Tau2 = TauAvg(ind2);
% Use linear interpolation to find current amplitude to produce 100Nm,
% rounded to the nearest integer value.
TauDes = 100;
sat = round((I2-I1)/(Tau2-Tau1)*(TauDes-Tau1)+I1);

% Clear any unneccessary variables.
clear answer minstim maxstim stim time mask_enables message_type mname format_string
clear LCstd LCm thresh_time thresh_ind LCmax sat_time sat_ind LC clear h
clear ind1 ind2 I1 I2 Tau1 Tau2 TauDes
clc

%% Procedure 6) Isometric Contractions
% The saturation (sat) is a required input variable to this experimental
% procedure.
Ready(6)

% Start the passive torque and position arrays with the
% hyperflexion/hyperextension data.
TauPass = [TauHF TauHE];
qPass = [qHF qHE];

choice = 'Yes';
i = 3;
while strcmp(choice,'Yes')
    
    %%%%%%%%%%%%%%%%%
    % Build the Simulink model
    if i == 3
        mname = 'Simulink Models/IsometricContractions';
        % Open the Simulink model
        open_system(mname)
        questdlg('Are you prepared to begin the isometric contraction test?', ...
            'Ready?', ...
            'Yes','No','Yes');
        qc_build_model
    else
        answer = questdlg('Are you prepared to begin the isometric contraction test, or would you like to take a 2 minute rest period?', ...
            'Ready?', ...
            'Ready','Rest','Ready');
        if strcmpi(answer,'Rest')
            pause(120)
        end
    end
    % Run the Simulink model
    qc_connect_model
    pause(3)
    qc_start_model
    % Pause and wait for the model to finish
    pause(9)
    %%%%%%%%%%%%%%%%%
    
    % Store the results into the results array.
    % Store the resulting passive data.
    TauPass(i) = TauPP;
    qPass(i) = qPP;
    % Store the resulting isometric contraction data.
    % The measured passive torque is subtracted so that only the torque due
    % to muscle stimulation is measured.
    TauIsoCon(i-2) = TauIso-TauPP;
    qIsoCon(i-2) = qIso;
    
    evalc(['TauIsoCon' num2str(i-2) ' = TauAct;']);
    evalc(['uIsoCon' num2str(i-2) ' = uAct;']);
    
    choice = questdlg('Would you like to perform another isometric contraction test?', ...
        'Isometric Contractions', ...
        'Yes','No','Yes');
    
    i = i+1;
end
% Close the Simulink model
save_system
close_system

% Clear unnecessary variables.
clear TauPP qPP TauIso qIso i mname answer
clear ans choice Current_filt format_string mask_enables message_type
clc

% The outputs of this procedure are the following arrays:
% TauPass: Passive joint torque measurements.
% qPass: Joint angles of the passive joint torque measurements.
% TauIsoCon: Joint torque measurements during the isometric contractions
% qIsoCon: Joint angles of isometric contraction tests.

% Each isometric contraction test also outputs the load cell data for the
% whole test duration (TauAct). This data can be used to compute the
% activation time constant.

%% Procedure 7) Sinusoidal Stimulation
% The subject is given a known smooth stimlation that varies between the
% saturation and the threshold. The input signal and the encoder signals
% will be saved to be used in the parameter ID process.
Ready(7)

U1 = 0.15;
U2 = 0.5;

F1 = 8;
F2 = F1;
F3 = F1;

% Setup the models time parameters
tend = 4*F1;

%%%%%%%%%%%%%%%%%
mname = 'Simulink Models/SinusoidalStimulation';
% Open the Simulink model
open_system(mname)
% Build the Simulink model
qc_build_model
% Run the Simulink model
qc_connect_model
pause(3)
qc_start_model
% Pause and wait for the model to finish
pause(tend+5)
% Close the Simulink model
save_system
close_system
%%%%%%%%%%%%%%%%%

% Clear unnecessary variables.
clear mname Current_filt mask_enables format_string message_type
clear F1 F2 F3 U1 U2 tend
clc

% The outputs of the model are normalized stimulation, knee angle, and
% time. This data will be used to determine the force velocity parameter.

% These variables are called u, phi, and time in the workspace.

%% Save Collected Data
% Organize the data that has been collected in the above procedures and
% save it for processing.

stamp = datestr(clock,'yyyy_mm_dd_HH_MM');
fname = ['Results_' stamp];

load SData

% Save the file name to the subject IDs
evalc(['aa = isfield(SData.' SID ',''allres'');']);
if aa
    evalc(['NF = numel(SData.' SID '.allres);']);
else
    NF = 0;
end

evalc(['SData.' SID '.allres{NF+1} = fname;']);

% Update SData with the calculated thresh and sat values
evalc(['SData.' SID '.params' Leg(1) '.thresh = thresh;']);
evalc(['SData.' SID '.params' Leg(1) '.sat = sat;']);
save('SData','SData')

clear SData aa

% Save the resulting data
save(['Results\' fname])
 
 

 