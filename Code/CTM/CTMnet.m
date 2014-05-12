% CTMnet: returns a vector of densities for a link for the next time step
% using the CTM
%
% INPUTS
% xnlink: vector of densities at the current time step for current link
% clink: a link object (not mapLinks)
% DT: time step (in seconds)
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

function xn1link=CTMnet(xnlink,clink,DT,isApp,err_S,err_Q,err_R,Q_S,Q_Q,Q_R)

% Link parameters
noLanes=clink.noLanes;
rhoj=clink.rhoj;
qmax=clink.qmax;
vmax=clink.vmax;
trueqin=clink.trueInflow; % already accounts for numLanes
trueqout=clink.trueOutflow; % already accounts for numLanes
dx=clink.dx;

% Compute dt/dx, a constant for the link
C=(DT/3600)/dx;

% Initialize
xn1link=zeros(1,clink.noCells);
icflows=zeros(clink.noCells+1,1);

% Get intercell flows
for i=1:length(icflows)
    if i==1
        icflows(i)=trueqin;
    elseif i==length(icflows)
        icflows(i)=trueqout;
    else
        icflows(i)=riemann_q(xnlink(i-1),xnlink(i),rhoj,qmax,vmax,noLanes,...
            isApp,err_S,err_Q,err_R,Q_S,Q_Q,Q_R);
    end
end

% Compute new densities
for i=1:clink.noCells
    qin=icflows(i);
    qout=icflows(i+1);
    xn1link(i)=xnlink(i)+C*(qin-qout);
end
        

