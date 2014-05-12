% getPex_Diag: takes in the PGA and the curves and returns the corresponding
% probability of exceedance for the damage state for each bridge (DIAGNOSTICS)
%
% INPUTS
% PGA: vector of PGA values
% p: array containing fragility curves
% SF: scale factor
% PGAVec: vector of PGA values for fragility model (defines range)

function Pex=getPex_Diag(PGA,p,SF,PGAVec)

% Initialize
num=length(PGA);
Pex=zeros(length(PGA),3);

for b=1:num
    
    % PGA rounded to the precision of SaVec
    roundedSa=round(PGA(b)*SF)/SF;
    
    % Get the corresponding index. 
    ind=find(abs(PGAVec-roundedSa)<.0001);
    
    % Get the Pex values
    for i=1:3
        
        Pex(b,i)=p(i,ind);
        
    end
    
end