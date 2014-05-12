% attachSensors: attach the sensors onto their appropriate links and
% overwrite previous map obj
%
% INPUTS
% mapLinks: a mapLinks object
% mapSensors: a mapSensors object
% numLinks: integer for the number of sensors
% keySensors: a cell array of sensor IDs

function mapLinks=attachSensors(mapLinks,mapSensors,numSensors,keySensors)

for i=1:numSensors
    
    % Current key
    ckey=keySensors{i};
    
    % Current sensor
    csensor=mapSensors(ckey);
    
    % Get the key of the appropriate link
    ckeyL=csensor.link;
    
    % New link object
    nlink=mapLinks(ckeyL);
    
    if length(nlink.sensorObjs)>=1
        
        ind=length(nlink.sensorObjs)+1;
        nlink.sensorObjs(ind)=csensor;
        
    else
        
        nlink.sensorObjs=csensor;
        
    end
    
    % Overwrite
    mapLinks(ckeyL)=nlink;
end