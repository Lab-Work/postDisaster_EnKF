% forwardPredictEns: runs the CTM forward one step and returns a vector of
% densities at the next time step (stochastic model)
%
% INPUTS
% DT: time step (in seconds)
% xipk: vector of densities at the current time step (posterior)
% mapLinks: a mapLinks object
% mapNodes: a mapNodes object
% totCells: integer for the number of cells in the system
% leftBC: upstream ghost cell density
% rightBC: downstream ghost cell density
% isApp: integer that says whether or not this is an approximation (open
% loop or filter) (N=0,Y=1)
% err_S: mean error in sending function
% err_Q: mean error in max flow region
% err_R: mean error in receiving function
% Q_S: variance in sending function
% Q_Q: variance in max flow region
% Q_R: variance in receiving function
% noiseBC: vector containing the uncertainty on the boundary conditions

function ximk1=forwardPredictEns(DT,xipk,mapLinks,mapNodes,totCells,...
    leftBC,rightBC,isApp,err_S,err_Q,err_R,Q_S,Q_Q,Q_R,noiseBC)

numLinks=length(mapLinks);
numNodes=length(mapNodes);

% Update the nodes
mapLinks=updateNodes(mapNodes,numNodes,mapLinks,numLinks,xipk,isApp,...
    err_S,err_Q,err_R,Q_S,Q_Q,Q_R);

% Update the links
postSt=updateLinks(xipk,mapLinks,numLinks,totCells,DT,leftBC,rightBC,isApp,...
    err_S,err_Q,err_R,Q_S,Q_Q,Q_R,noiseBC);

% Transpose
ximk1=postSt';
    

