% updateQ: creates new link objects and their associated keys
%
% INPUTS
% mapLinks: a mapLinks object
% cnode: a node object (not mapNodes)
% x: vector of densities at teh current time step
% isApp: integer that says whether or not this is an approximation (open
% loop or filter) (N=0,Y=1)
%
% The following are inputs but not used in the deterministic model
% 
% err_S: mean error in sending function
% err_Q: mean error in max flow region
% err_R: mean error in receiving function
% Q_S: variance in sending function
% Q_Q: variance in max flow region
% Q_R: variance in receiving function

function [nLinks,lkeys]=updateQ(mapLinks,cnode,x,isApp,...
    err_S,err_Q,err_R,Q_S,Q_Q,Q_R)

% For the incoming and outgoing links of the node, determine the keys, link
% objects, and indices of the critical cells

% Incoming links
key_in=cnode.links_in;
lin=mapLinks(key_in);
lin_ind=lin.endCell;

% Outgoing links
key_out=cnode.links_out;
lout=mapLinks(key_out);
lout_ind=lout.startCell;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Getting/assigning q's

rhoj1=lin.rhoj;
qmax1=lin.qmax;
vmax1=lin.vmax;

rhoj2=lout.rhoj;
qmax2=lout.qmax;
vmax2=lout.vmax;

% Getting number of lanes
lanesIn=lin.noLanes;
lanesOut=lout.noLanes;

% Compute sending and receiving functions
S1=sending_q(x(lin_ind),rhoj1,qmax1,vmax1,lanesIn,err_S,err_Q,Q_S,Q_Q,isApp);
R2=receiving_q(x(lout_ind),rhoj2,qmax2,vmax2,lanesOut,err_R,err_Q,Q_R,Q_Q,isApp);

% Compute flow
q12=min(S1,R2);

% Assign to link
lout.trueInflow=q12;
lin.trueOutflow=q12;

nLinks=[lout lin];
lkeys={key_out;key_in};