% updateMergeQ: creates new link objects and their associated keys
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

function [nLinks,lkeys]=updateMergeQ(mapLinks,cnode,x,isApp,...
    err_S,err_Q,err_R,Q_S,Q_Q,Q_R)

% For the incoming and outgoing links of the node, determine the keys, link
% objects, and indices of the critical cells

% Outgoing links
key_out=cnode.links_out;
lout=mapLinks(key_out);
lout_ind=lout.startCell;

% Incoming links
str_in=cnode.links_in;
key_in=regexp(str_in,'\W+','split');

lin=link.empty(0,length(key_in));
lin_ind=zeros(1,length(key_in));

for i=1:length(key_in)
    
    lin(i)=mapLinks(key_in{i});
    lin_ind(i)=lin(i).endCell;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate beta values based on number of incoming lanes. This beta value
% essentially allocates how much inflow is coming from each link into the
% node.

totInLanes=0;

for i=1:length(lin)
    
    totInLanes=totInLanes+lin(i).noLanes;
    
end

for i=1:length(lin)
    
    lin(i).beta=lin(i).noLanes/totInLanes;
    
end

% Getting/assigning q's

rhoj1=lin(1).rhoj;
qmax1=lin(1).qmax;
vmax1=lin(1).vmax;
rhoj2=lin(2).rhoj;
qmax2=lin(2).qmax;
vmax2=lin(2).vmax;

rhoj3=lout.rhoj;
qmax3=lout.qmax;
vmax3=lout.vmax;

% Getting number of lanes
lanesIn1=lin(1).noLanes;
lanesIn2=lin(2).noLanes;
lanesOut=lout.noLanes;

% Compute sending and receiving functions
S1=sending_q(x(lin_ind(1)),rhoj1,qmax1,vmax1,lanesIn1,err_S,err_Q,Q_S,Q_Q,isApp);
S2=sending_q(x(lin_ind(2)),rhoj2,qmax2,vmax2,lanesIn2,err_S,err_Q,Q_S,Q_Q,isApp);
R3=receiving_q(x(lout_ind),rhoj3,qmax3,vmax3,lanesOut,err_R,err_Q,Q_R,Q_Q,isApp);

% Compute inflows
if S1+S2>R3
    q13=min(S1,lin(1).beta*R3);
    q23=min(S2,lin(2).beta*R3);
    
else
    q13=min(S1,R3);
    q23=min(S2,R3);
    
end

% Assign to links
lin(1).trueOutflow=q13;
lin(2).trueOutflow=q23;
lout.trueInflow=q13+q23;

nLinks=[lin(1) lin(2) lout];
lkeys=[key_in key_out];
