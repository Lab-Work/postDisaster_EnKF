% composeEnsToMat: Create the A matrix of forecasted ensemble members
% according to Evensen
%
% INPUTS
% xi: 4D array of density ensembles (prior or posterior)
% totCells: integer for number of cells in system
% noEn: number of ensembles
% timeInd: the desired time index

function A=composeEnsToMat(xi,totCells,noEn,timeInd)

% Initialize
A=zeros(totCells,noEn);

for i=1:noEn
    
    % Put ensembles into matrix
    A(:,i)=xi(:,1,i,timeInd);
    
end