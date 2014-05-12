% sending_q: sending function with qmax and rhoj as inputs for triangular
% fundamental diagram
%
% INPUTS:
% rho: density of current cell
% rhoj: jam density of link
% qmax: max flow of link
% vmax: max speed of cars in link
% numLanes: integer for number of lanes of link
% err_S: mean error in sending function
% err_Q: mean error in max flow region
% Q_S: variance in sending function
% Q_Q: variance in max flow region
% isApp: integer that says whether or not this is an approximation (open
% loop or filter) (N=0,Y=1)

function S=sending_q(rho,rhoj,qmax,vmax,numLanes,err_S,err_Q,Q_S,Q_Q,isApp)

% Account for number of lanes
qmax=qmax*numLanes;
rhoj=rhoj*numLanes;

rhoc=qmax/vmax;

% Deterministic
if rho>=rhoc
    S=qmax;
else
    S=vmax*rho;
end

% Noise (if applicable)
if isApp==1
    if rho>=rhoc
        err=normrnd(err_Q,sqrt(Q_Q));
    else
        err=normrnd(err_S,sqrt(Q_S));
    end
else
    err=0;
end

S=S+err;

S=max(S,0);
