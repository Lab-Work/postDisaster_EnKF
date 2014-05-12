% alterParams: Slightly perturbs the values of the lane parameters in the
% true model. It gives back a new valueLinks array and it gives the value
% of vmax drawn for each link in the network.
% 
% INPUTS:
% valueLinks: a valueLinks array

function [vmaxVec,nvalueLinks]=alterParams(valueLinks)

% Initialization
nvalueLinks=valueLinks;
numLinks=length(valueLinks);

% Noise (std dev %)
maxCapNoise=5;
vmaxNoise=3;
rhojNoise=5;

for i=1:numLinks
    
    % Get current link
    nlink=valueLinks{i};
    
    % Alter parameters
    nlink.maxCap=normrnd(nlink.maxCap,(maxCapNoise/100)*nlink.maxCap);
    nlink.vmax=normrnd(nlink.vmax,(vmaxNoise/100)*nlink.vmax);
    nlink.maxRhoj=normrnd(nlink.maxRhoj,(rhojNoise/100)*nlink.maxRhoj);
    
    % Replace link and populate vmax vector
    nvalueLinks{i}=nlink;
    vmaxVec(i)=nlink.vmax;
    
end

        
        

