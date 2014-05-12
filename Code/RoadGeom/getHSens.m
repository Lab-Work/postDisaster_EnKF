% getHSens: get the measurement matrix and dimensions
%
% INPUTS
% mapSensors: a mapSensors object
% keySensors: a cell array of sensor IDs
% totCells: number of cells in system
% numSensors: integer for number of sensors
% 
% OUTPUTS
% H: measurement matrix
% kk: number of rows of H
% n: number of columns of H

function [H kk n]=getHSens(mapSensors,keySensors,totCells,numSensors)

% Initialize
H=zeros(numSensors,totCells);

for i=1:numSensors
    
    % Current key
    ckey=keySensors{i};
    
    % Current sensor
    csensor=mapSensors(ckey);
    
    % Global position of sensor
    gPos=csensor.globalCell;
    
    % Update H
    H(i,gPos)=1;
  
end

% Get dimensions
kk=size(H,1);
n=size(H,2);