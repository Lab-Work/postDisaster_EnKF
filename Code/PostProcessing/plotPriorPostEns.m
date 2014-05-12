% plotPriorPostEns: plot the prior and posterior ensembles for all time
% steps showing true measurements and ensemble measurements
%
% INPUTS
% xi_m: prior ensemble array
% xi_p: posterior ensemble array
% D_array: cell array that contains D matrix at every time step
% sortInd: vector that contains integer indices of sensor locations
% y: array that contains measuerments
% start: start time step
% fin: end time step
% noEn: integer for the number of ensembles
% measTime: vector of times
% pDuration: pause duration (may or may not be used)
% 
% NOTE: There might be a time lag between time steps

function plotPriorPostEns(xi_m,xi_p,D_array,sortInd,y,start,fin,noEn,measTime,pDuration)

figure

for k=start:fin 
    
    clf
        
    % Prior
    subplot(1,2,1)
    hold on
    for i=1:noEn
        plot(xi_m(:,1,i,k),'g')
    end
    
    % Add measurements
    boxplot(D_array{k}','positions',sortInd,'widths',.75)
    title(['Prior, k=' num2str(k) ', T (s)=' num2str(measTime(k))...
        ', y(k)= ' num2str(y(k,:))]);
    xlabel('Cell')
    ylabel('Density')
    set(gca,'XTickLabelMode','auto','XTickMode','auto')
    grid on
    axis auto
    
    hold off
    
    % Posterior
    subplot(1,2,2)
    hold on
    for i=1:noEn
        plot(xi_p(:,1,i,k),'g')
    end
    
    % Add measurements
    boxplot(D_array{k}','positions',sortInd,'widths',.75)
    title(['Post, k=' num2str(k) ', T (s)=' num2str(measTime(k))...
        ', y(k)= ' num2str(y(k,:))]);
    xlabel('Cell')
    ylabel('Density')
    set(gca,'XTickLabelMode','auto','XTickMode','auto')
    grid on
    axis auto
    
    hold off

    pause(pDuration);
    
    % Manually operate pause
%     pause
    


end