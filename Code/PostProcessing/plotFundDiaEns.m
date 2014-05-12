% plotFundDiaEns: Plots the fundamental diagram for each ensemble in a
% desired range.
%
% INPUTS
% mapLinks_array: a mapLinks array
% keyLinks: a cell array of link IDs
% numLinks: integer for the number of links
% startEn: the starting ensemble
% finEn: the end ensemble
% pDuration: duration if the pause (if applicable)


function plotFundDiaEns(mapLinks_array,keyLinks,numLinks,startEn,finEn,pDuration)

% True parameters (hardcoded)
qmax_t=2000;
vmax_t=110;
rhoj_t=125;

figure

for i=startEn:finEn
    
    mapLinks=mapLinks_array{i};
    
    for j=1:numLinks
        
        % Current key
        ckey=keyLinks{j};
        
        % Current link
        clink=mapLinks(ckey);
        
        % NOISY FUNDAMENTAL DIAGRAM
        
        qmax=clink.qmax;
        vmax=clink.vmax;
        rhoj=clink.rhoj;
        
        x1=0;
        x2=qmax/vmax;
        x3=rhoj;
        
        y1=0;
        y2=qmax;
        y3=0;
        
        subplot(2,3,j)
        line([x1 x2 x3],[y1 y2 y3]);
        
        hold on
        
        % TRUE FUNDAMENTAL DIAGRAM
        
        x1_t=0;
        x2_t=qmax_t/vmax_t;
        x3_t=rhoj_t;
        
        y1_t=0;
        y2_t=qmax_t;
        y3_t=0;
        
        line([x1_t x2_t x3_t],[y1_t y2_t y3_t],'Color',[1 0 0]);
        
        title({['Link: ' ckey, ', Ens: ' num2str(i)],...
            ['v: ' num2str(vmax) ', q: ' num2str(qmax) ', \rho_j: ' num2str(rhoj)],...
            ['\rho_c: ' num2str(qmax/vmax) ', w: ' num2str(qmax/(rhoj-qmax/vmax))]})
        xlabel('\rho')
        ylabel('q')
        legend('Ensemble Draw','True')
        
    end
    
    pause(pDuration);
    
    % Manually control pause
%     pause

    clf
    
end