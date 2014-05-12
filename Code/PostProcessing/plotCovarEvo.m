% plotCovarEvo: plots the evolution of the covariance matrices for open
% loop, prior, and posterior solutions in contour form.
% 
% INPUTS
% Pk_p2: error covariance matrix for open loop sol
% Pk_m: error covariance matrix for filtered sol (prior)
% Pk_p: error covariance matrix for filtered sol (posterior)
% start: start time step
% fin: end time step
% measTime: vector of times
% pDuration: pause duration (may or may not be used)

function plotCovarEvo(Pk_p2,Pk_m,Pk_p,start,fin,measTime,pDuration)

figure

for i=start:fin
    
    % Open Loop
    subplot(2,2,1)
    imagesc(Pk_p2(:,:,i))
    title(['Open Loop, Step= ' num2str(measTime(i)/5+1) ', Time (s) = ' num2str(measTime(i))]);
    caxis([0 100]);
    colorbar;
    
    % Before Meas
    subplot(2,2,2)
    imagesc(Pk_m(:,:,i))
    title(['Prior, Step= ' num2str(measTime(i)/5+1) ', Time (s) = ' num2str(measTime(i))]);
    caxis([0 100]);
    colorbar;
    
    % After Meas
    subplot(2,2,3)
    imagesc(Pk_p(:,:,i))
    title(['Posterior, Step= ' num2str(measTime(i)/5+1) ', Time (s) = ' num2str(measTime(i))]);
    caxis([0 100]);
    colorbar;
    
    pause(pDuration);
    
    % Manually control pause
%     pause
    
end