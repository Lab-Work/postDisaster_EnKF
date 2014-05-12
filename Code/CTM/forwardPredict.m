% forwardPredict: runs the CTM forward one step and returns a vector of
% densities at the next time step (deterministic model)
%
% INPUTS
% DT: time step (in seconds)
% xk: vector of densities at teh current time step
% mapLinks: a mapLinks object
% mapNodes: a mapNodes object
% totCells: integer for the number of cells in the system
% leftBC: upstream ghost cell density
% rightBC: downstream ghost cell density
% isFilter: integer that says whether or not this is the estimator
% (N=0,Y=1)
% noiseBC: vector containing the uncertainty on the boundary conditions
%
% The following are inputs but not used in the deterministic model
% 
% err_S: mean error in sending function
% err_Q: mean error in max flow region
% err_R: mean error in receiving function
% Q_S: variance in sending function
% Q_Q: variance in max flow region
% Q_R: variance in receiving function

function xk1=forwardPredict(DT,xk,mapLinks,mapNodes,totCells,leftBC,rightBC,...
    isFilter,err_S,err_Q,err_R,Q_S,Q_Q,Q_R,noiseBC)

numLinks=length(mapLinks);
numNodes=length(mapNodes);

% Get the true outflows and true inflows from each node and assign these as
% inflows and outflows on the link objects
mapLinks=updateNodes(mapNodes,numNodes,mapLinks,numLinks,xk,isFilter,...
    err_S,err_Q,err_R,Q_S,Q_Q,Q_R);

% Update the densities using the CTM
xk1=updateLinks(xk,mapLinks,numLinks,totCells,DT,leftBC,rightBC,isFilter,...
    err_S,err_Q,err_R,Q_S,Q_Q,Q_R,noiseBC);
