% getPex: takes in the PGA and the curves and returns the corresponding
% probability of exceedance for the damage state for each bridge
%
% INPUTS
% PGA: vector of PGA values
% p: array containing fragility curves
% SF: scale factor
% PGAVec: vector of PGA values for fragility model (defines range)

function Pex=getPex(PGA,p,SF,PGAVec)

% Initialize
numBridges=length(PGA);
numCurves=3; % hardcoded
Pex=zeros(length(PGA),numCurves);

for i=1:numBridges
    
    % PGA rounded to the precision of SaVec
    roundedSa=round(PGA(i)*SF)/SF;
    
    % Get the corresponding index. 
    ind=find(abs(PGAVec-roundedSa)<.0001);
    
    % Get the Pex values
    for j=1:3
        
        Pex(i,j)=p(j,ind,i);
        
    end
    
end