% getSensorPos: returns sensors global position in cell network in new map
% object
%
% INPUTS
% mapSensors: a mapSensors object
% mapLinks: a mapLinks object
% numSensors: integer for the number of sensors
%
% OUTPUTS
% mapSensors: modified mapSensors object
% orderedKeySens: a keySensors map obj whose keys are ordered by where the
% sensors are in the system
% sortInd: a vector that contains integer cell indices of the sensor
% locations sorted numerically

function [mapSensors,orderedKeySens,sortInd]=getSensorPos(mapSensors,...
    mapLinks,numSensors)

% Initialize
gvec=zeros(1,numSensors);
orderedKeySens=cell(1,numSensors);

% Get keys
keySensors=keys(mapSensors);

for i=1:numSensors
    
     % Current key
    ckey=keySensors{i};
    
    % Current sensor
    csensor=mapSensors(ckey);
    
    % Link associated to current sensor
    linkName=csensor.link;
    
    % Call link object
    clink=mapLinks(linkName);
    
    % Get necessary values
    off=csensor.offset;
    startC=clink.startCell;
    dx=clink.dx;
    
    % Compute global position of the current sensor obj
    if off==0
        localCell=1;
    else
        localCell=ceil(off/dx);
    end
    
    globalPos=startC+localCell-1;
    
    csensor.globalCell=globalPos;
    
    % Overwrite old sensor obj
    mapSensors(ckey)=csensor;
    
    gvec(i)=globalPos;
    
end

% Rearrange

[sortInd,order]=sort(gvec);

for i=1:numSensors
    
    orderedKeySens(i)=keySensors(order(i));
    
end