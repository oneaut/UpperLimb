%% Data Corrections
% This script will perform any necessary corrections on the data that was
% collected during the parameter estimation procedure.

switch fname
    
    case 'Results_2016_04_20_10_32'         % Participant C1
        % Pushed when I should have pulled during the push/pull test
        tmp = TauHE;
        TauHE = TauHF;
        TauHF = tmp;
        tmp = qHE;
        qHE = qHF;
        qHF = tmp;
        
        % Sign of all of the load cell data is incorrect
        LCsgn = -1;
        
    case 'Results_2016_04_20_17_56'         % Participant C1
        % Sign of all of the load cell data is incorrect
        LCsgn = -1;
        
    case 'Results_2016_05_12_14_19'         % Participant C1
        % Sign of all of the load cell data is incorrect
        LCsgn = -1;
        
    case 'Results_2016_05_12_14_45'         % Participant C1
        % Sign of all of the load cell data is incorrect
        LCsgn = -1;
        
    case 'Results_2016_04_26_15_34'         % Participant C2
        % Sign of all of the load cell data is incorrect
        LCsgn = -1;
        
    case 'Results_2016_04_26_16_18'         % Participant C2
        % Sign of all of the load cell data is incorrect
        LCsgn = -1;
        
    case 'Results_2016_05_19_11_38'         % Participant C2
        % Sign of all of the load cell data is incorrect
        LCsgn = -1;
        % Leg was not lifted early enough during pendulum test
        qPend = qPend(2001:end);
        TauPend = TauPend(2001:end);
        % Load cell was faulty during pendulum test
        TauPend = zeros(size(TauPend));
        % The load cell malfunctioned during the hyper extension/flexion
        % tests. These correction data points are pulled from the tests on
        % 5/20/2016 at 12:00PM.
        TauPass(1) = 20.41;
        TauPass(2) = -13.6;
        
    case 'Results_2016_05_19_12_00'         % Participant C2
        % Sign of all of the load cell data is incorrect
        LCsgn = -1;
        
end
