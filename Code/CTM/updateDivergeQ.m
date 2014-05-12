% updateDivergeQ: creates new link objects and their associated keys
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

function [nLinks,lkeys]=updateDivergeQ(mapLinks,cnode,x,isApp,...
    err_S,err_Q,err_R,Q_S,Q_Q,Q_R)

% For the incoming and outgoing links of the node, determine the keys, link
% objects, and indices of the critical cells

% Incoming links
key_in=cnode.links_in;
lin=mapLinks(key_in);
lin_ind=lin.endCell;

% Outgoing links
str_out=cnode.links_out;
key_out=regexp(str_out,'\W+','split');

lout=link.empty(0,length(key_out));
lout_ind=zeros(1,length(key_out));

for i=1:length(key_out)
    
    lout(i)=mapLinks(key_out{i});
    lout_ind(i)=lout(i).startCell;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate alpha values based on number of incoming lanes. This alphs value
% essentially allocates how much outflow is coming from each link out of the
% node.

totOutLanes=0;

for i=1:length(lout)
    
    totOutLanes=totOutLanes+lout(i).noLanes;
    
end

for i=1:length(lout)
    
    lout(i).alpha=lout(i).noLanes/totOutLanes;
    
end

% Getting/assigning q's

rhoj1=lin.rhoj;
qmax1=lin.qmax;
vmax1=lin.vmax;

rhoj2=lout(1).rhoj;
qmax2=lout(1).qmax;
vmax2=lout(1).vmax;
rhoj3=lout(2).rhoj;
qmax3=lout(2).qmax;
vmax3=lout(2).vmax;

% Getting number of lanes
lanesIn=lin.noLanes;
lanesOut1=lout(1).noLanes;
lanesOut2=lout(2).noLanes;

% Compute sending and receiving functions
S1=sending_q(x(lin_ind),rhoj1,qmax1,vmax1,lanesIn,err_S,err_Q,Q_S,Q_Q,isApp);
R2=receiving_q(x(lout_ind(1)),rhoj2,qmax2,vmax2,lanesOut1,err_R,err_Q,Q_R,Q_Q,isApp);
R3=receiving_q(x(lout_ind(2)),rhoj3,qmax3,vmax3,lanesOut2,err_R,err_Q,Q_R,Q_Q,isApp);

% Compute outflows
q12=min(lout(1).alpha*S1,R2);
q13=min(lout(2).alpha*S1,R3);

% Assign to links
lout(1).trueInflow=q12;
lout(2).trueInflow=q13;
lin.trueOutflow=q12+q13;

nLinks=[lout(1) lout(2) lin];
lkeys=[key_out key_in];