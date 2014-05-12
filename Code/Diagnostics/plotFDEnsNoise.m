% plotFDEnsNoise: plots true fundamental diagram and sample points drawn
% from a distribution.

clear all;close all;clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TRUE PARAMS
t_qmax=2000; % veh/h/lane
t_vmax=110; % km/h
t_rj=125; %veh/km

% Noise percentages
qNoise=5;
vNoise=3;
rjNoise=5;

noEn=100; % number of ens
noPtperEns=200; % noisy samples are ens

% Type of noise ('LN' or 'N')
typeNoise='N';

% Sending and receiving function uncertainty
mu_s=150;
mu_q=200;
mu_r=300;
var_s=100^2;
var_q=150^2;
var_r=200^2;

% Lognormal params
err_S=log(mu_s^2/sqrt(mu_s^2+var_s)); % lambda
err_Q=log(mu_q^2/sqrt(mu_q^2+var_q)); % lambda
err_R=log(mu_r^2/sqrt(mu_r^2+var_r)); % lambda
Q_S=log(1+var_s/mu_s^2); % xi^2
Q_Q=log(1+var_q/mu_q^2); % xi^2
Q_R=log(1+var_r/mu_r^2); % xi^2

% Normal noise params (standard deviation)
s_norm=50;
r_norm=150;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

hold on

for i=1:noEn
    
    q=normrnd(t_qmax,(qNoise/100)*t_qmax);
    v=normrnd(t_vmax,(vNoise/100)*t_vmax);
    rj=normrnd(t_rj,(rjNoise/100)*t_rj);

    % Compute w and rc
    w=q/(rj-q/v);
    rc=q/v;
    
    x1=0;
    x2=rc;
    x3=rj;
    
    y1=0;
    y2=q;
    y3=0;
    
    line([x1 x2 x3],[y1 y2 y3]);
    
    hold on
    
    % NOISY FUNDAMENTAL DIAGRAM
    
    a=0;
    b=rj;
    
    randrho=a+(b-a).*rand(noPtperEns,1);
    val=zeros(noPtperEns,1);
    fund_noise=zeros(noPtperEns,1);
  
    % Obtain indices of freeflow and congested vectors
    ffID=randrho<rc;
    cID=randrho>=rc;
    
    % Draw random vectors
    val(ffID)=randrho(ffID)*v;
    val(cID)=w*(rj-randrho(cID));   
    if isequal(typeNoise,'LN')
      % do the LN noise
      errFF=-1*lognrnd(err_S,sqrt(Q_S),sum(ffID),1);
      errC=-1*lognrnd(err_R,sqrt(Q_R),sum(cID),1);
    else
      % do the other one
      errFF=normrnd(0,s_norm,sum(ffID),1);
      errC=normrnd(0,r_norm,sum(cID),1);
    end  
    fund_noise(ffID)=val(ffID)+errFF;
    fund_noise(cID)=val(cID)+errC;
           
    plot(randrho,fund_noise,'k.','MarkerSize',1);
    
    hold on
    
    % TRUE FUNDAMENTAL DIAGRAM

    t_rhoc=t_qmax/t_vmax;
    t_rhoj=t_rj;
    
    x1_t=0;
    x2_t=t_rhoc;
    x3_t=t_rhoj;
    
    y1_t=0;
    y2_t=t_qmax;
    y3_t=0;
    
    line([x1_t x2_t x3_t],[y1_t y2_t y3_t],'LineWidth',2,'Color',[1 0 0]);
    
    title({['No ens: ' num2str(noEn) ', No Pts per Ens: ' num2str(noPtperEns)]...
        ['Fund noises (pct): qmax: ' num2str(qNoise) ', vmax: ' num2str(vNoise) ', rj: ' num2str(rjNoise)]});
    xlabel('\rho')
    ylabel('q(\rho)')
    ylim([0 3000])
    legend('Ensemble Draw','Noise','True')
    
    % Comment to plot all ensembles at once
%     pause
%     clf
    
end