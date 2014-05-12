% getEnsMeasNoiseEQ: takes the measurement, computes an ensemble noise, and
% returns both the ensemble noise vector and the noisy measurement ensemble
% vector when Heq (modified H matrix) is used.
%
% INPUTS
% meas: the measurement vector at the desired times step
% mapSensors: a mapSensors object
% keySensors: a cell array of sensor IDs
% sensInd: a vector of indices corresponding to sensors in Heq

function [ensMeas ensMeasNoise]=getEnsMeasNoiseEQ(meas,mapSensors,keySensors,sensInd)

% Initialize
ensMeasNoise=zeros(length(meas),1);

for i=1:length(sensInd)
    
    % Current key
    ckey=keySensors{sensInd(i)};
    
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
