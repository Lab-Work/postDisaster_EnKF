% getEnsMeasNoise: takes the measurement, computes an ensemble noise, and
% returns both the ensemble noise vector and the noisy measurement ensemble
% vector
%
% INPUTS
% meas: the measurement vector at the desired times step
% mapSensors: a mapSensors object
% keySensors: a cell array of sensor IDs

function [ensMeas ensMeasNoise]=getEnsMeasNoise(meas,mapSensors,keySensors)

% Initialize
ensMeasNoise=zeros(length(meas),1);

for i=1:length(meas)
    
    % Current key
    ckey=keySensors{i};
    
    % Current link
    csens=mapSensors(ckey);
    
    % Get necessary params
    sensMu=csens.mu;
    sensSig=csens.sigma;
   
    % Get noise vector
    ensMeasNoise(i)=normrnd(sensMu,sensSig);
    
end

% Noisy observation
ensMeas=meas+ensMeasNoise;
