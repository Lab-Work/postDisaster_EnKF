% getProbVec_Diag: takes into account the probability of exceeedances and
% outputs the probability of damage states as a vector (or matrix for
% multiple bridges) (DIAGNOSTICS)
%
% INPUTS
% Pex: probability of exceedance matrix
%
% NOTE: This is currently hardcoded for 4 damage states

function ProbDS=getProbVec_Diag(Pex)

% Initialize
num=size(Pex,1);
ProbDS=zeros(num,4);

for i=1:num
    
    ProbDS(i,1)=Pex(i,3);
    ProbDS(i,2)=Pex(i,2)-Pex(i,3);
    ProbDS(i,3)=Pex(i,1)-Pex(i,2);
    ProbDS(i,4)=1-Pex(i,1);
    
end