% riemann_q: computes the minimum of sending and receiving functions (flow)
% across an interface given a link
%
% INPUTS
% rho1: density of upstream cell
% rho2: density of downstream cell
% rhoj: jam density of link
% qmax: max flow of link
% vmax: max speed of cars in link
% numLanes: integer for number of lanes of link
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

function q=riemann_q(rho1,rho2,rhoj,qmax,vmax,numLanes...
    ,isApp,err_S,err_Q,err_R,Q_S,Q_Q,Q_R)

% Account for number of lanes
qmax=qmax*numLanes;
rhoj=rhoj*numLanes;

rhoc=qmax/vmax;

% Sending
s=(rho1>=rhoc).*(qmax)+...
    (rho1<rhoc).*(vmax*rho1);

% Receiving
r=(rho2>=rhoc).*((qmax/(rhoj-rhoc)).*(rhoj-rho2))+...
    (rho2<rhoc).*(qmax);

% Draw error
if isApp==1    
    if rho1<rhoc
        s_Err=normrnd(err_S,sqrt(Q_S));
    else
        s_Err=normrnd(err_Q,sqrt(Q_Q));

    end
    
    if rho2<rhoc
        r_Err=normrnd(err_Q,sqrt(Q_Q));
    else
        r_Err=normrnd(err_R,sqrt(Q_R));
    end
    
else
    s_Err=0;
    r_Err=0;
end

% Exclude negative values of flow
s=s+s_Err;
r=r+r_Err;
q=max(min(s,r),0);