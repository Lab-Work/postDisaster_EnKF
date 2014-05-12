% getHeq: get the modified measurement matrix and the indices of the
% sensors
%
% INPUTS
% mapSensors: a mapSensors object
% keySensors: a cell array of sensor IDs
% totCells: number of cells in system
% numSensors: integer for number of sensors
% pf_sens: the probability of failure of a single sensor
% 
% OUTPUTS
% Heq: modified measurement matrix
% sensInd: indices of the sensors that havent failed

function [Heq sensInd]=getHeq(mapSensors,keySensors,totCells,numSensors,pf_sens)

% Initialize
Heq=zeros(1,totCells);

% Keep track of number of sensors
ind=1;

for i=1:numSensors
    
    draw=rand;
    
    % Doesnt fail
    if draw>pf_sens
        
        % Current key
        ckey=keySensors{i};
        
        % Current sensor
        csensor=mapSensors(ckey);
        
        % Global position of sensor
        gPos=csensor.globalCell;
        
        % Update H
        Heq(ind,gPos)=1;
        
        % Populate and increment index
        sensInd(ind)=i;
        ind=ind+1;
        
    end
    
end