% inpSensors: Creates traffic sensors by giving them an ID. Also associated
% with the sensors is the link it is on, the type of sensor it is ('dens'),
% how far downstream from the beginning of the link it is, and the accuracy
% of the sensor.
%
% INPUTS
% numSensors: integer that tells the number of sensors in the network

function [valueSensors,keySensors]=inpSensors(numSensors)

% Initialize
valueSensors=cell(1,numSensors);
keySensors=cell(1,numSensors);

prompt={'ID','Associated Link','Sensor Type (dens,etc)','Offset',...
    'mu','sigma'};

for i=1:numSensors
    
    inp=inputdlg(prompt,['Sensor ' num2str(i)],1,{['s',num2str(i)],'0','dens','0',...
        '0','5'});
    
    s=sensor;
    
    s.ID=inp{1};
    s.link=inp{2};
    s.type=inp{3};
    s.offset=str2num(inp{4});
    s.mu=str2num(inp{5});
    s.sigma=str2num(inp{6});
    
    keySensors{i}=s.ID;
    valueSensors{i}=s;
    
end