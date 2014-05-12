% receiving_q_neg: receiving function with qmax and rhoj as inputs for
% triangular fundamental diagram
%
% INPUTS:
% rho: density of current cell
% rhoj: jam density of link
% qmax: max flow of link
% vmax: max speed of cars in link
% numLanes: integer for number of lanes of link
% err_R: mean error in receiving function
% err_Q: mean error in max flow region
% Q_R: variance in receiving function
% Q_Q: variance in max flow region
% isApp: integer that says whether or not this is an approximation (open
% loop or filter) (N=0,Y=1)
%
% NOTE: This is a special receiving function which allows for negative values

function R=receiving_q_neg(rho,rhoj,qmax,vmax,numLanes,err_R,err_Q,Q_R,Q_Q,isApp)

% Account for number of lanes
qmax=qmax*numLanes;
rhoj=rhoj*numLanes;

rhoc=qmax/vmax;

% Deterministic
if rho<=rhoc
    R=qmax;
else
    % Special code if rhoj=0. Otherwise you will get NaN
    if rhoj==0
        R=0;
    else
        R=(qmax/(rhoj-rhoc))*(rhoj-rho);
    end
end

% Noise (if applicable)
if isApp==1
    if rho>=rhoc
        err=normrnd(err_R,sqrt(Q_R));
    else
        err=normrnd(err_Q,sqrt(Q_Q));
    end
else
    err=0;
end

R=R+err;