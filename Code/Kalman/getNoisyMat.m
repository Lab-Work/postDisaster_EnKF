% getNoisyMat: creates a noise matrix that is added to the true solution to
% get synthetic measurements
%
% INPUTS
% mapSensors: a mapSensors object
% numSensors: integer for the number of sensors in the system
% keySensors: a cell array of sensor IDs
% measTime: vector of time steps (in seconds)

function noisyMat=getNoisyMat(mapSensors,numSensors,keySensors,measTime)

% Initialize
noisyMat=zeros(length(measTime),numSensors);

for i=1:numSensors
    
    % Get the key
    ckey=keySensors{i};
    
    % Get the sensor
    csensor=mapSensors(ckey);
    
    % Extract relevant parameters from sensor object
    cmu=csensor.mu;
    csig=csensor.sigma;
    
    % Create noisy vector
    noisyMat(:,i)=normrnd(cmu,csig,length(measTime),1);
    
end