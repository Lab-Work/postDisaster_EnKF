% getyMeas: computes synthetic measurements
%
% INPUTS
% H: measurement matrix
% xMatNoise: matrix of true solution densities (row: time, col: cell)
% measTime: vector of time steps (in seconds)
% numSensors: integer for the number of sensors in the system
% noisyMat: perturbation matrix of true solution

function y=getyMeas(H,xMatNoise,measTime,numSensors,noisyMat)

% Initialize
y=zeros(length(measTime),numSensors);

for i=1:length(measTime)
    
    % True solution
    x=xMatNoise(i,:)';
    
    % Perturbation
    v=noisyMat(i,:)';
    
    % Synthetic measurement
    y(i,:)=H*x+v;
    
end