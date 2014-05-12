% computeBEEQ: computes the BEEQ value for a given estimator
% 
% INPUTS
% xMatTrue: true solution
% xMatOpenLp: open loop solution
% xMatAnalyzed: estiamted solution

function BEEQ=computeBEEQ(xMatTrue,xMatOpenLp,xMatAnalyzed)

K=size(xMatTrue,1);
numCells=size(xMatTrue,2);

% Reshape the vectors in such a way to omit the boundary cells
trueSol=reshape(xMatTrue(:,1:numCells-2),1,K*(numCells-2));
openSol=reshape(xMatOpenLp(:,1:numCells-2),1,K*(numCells-2));
analyzedSol=reshape(xMatAnalyzed(:,1:numCells-2),1,K*(numCells-2));

% Compute BEEQ
numDiff=analyzedSol-trueSol;

num=norm(numDiff);

denomDiff=openSol-trueSol;

denom=norm(denomDiff);

BEEQ=num/denom;