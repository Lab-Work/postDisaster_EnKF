% plotCompCovar: comparison plot between true, est, and data solutions
% (posterior, matrix form). Also plots covariance bands with standard dev
% of a given value
% 
% INPUTS
% xMatNoise: true solution 
% xMatp: filtered solution
% xMatp2: open loop solution
% Pk_p: error covariance matrix for filtered sol
% Pk_sm: error covariance matrix for smoothed sol
% start: start time step
% fin: end time step
% totCells: number of cells in system
% measTime: vector of times
% sd: number of standard deviations of band
% pDuration: pause duration (may or may not be used)

function plotCompCovar(xMatNoise,xMatp,xMatp2,Pk_p,start,fin,totCells,...
    measTime,sd,pDuration)

% Reshape the matrices in such a way to omit the boundary cells
xMatNoise=xMatNoise(:,1:totCells-2);
xMatp2=xMatp2(:,1:totCells-2);
xMatp=xMatp(:,1:totCells-2);
Pk_p=Pk_p(1:totCells-2,1:totCells-2,:);

% Organize into a master array
matMaster(:,:,1)=xMatNoise;
matMaster(:,:,2)=xMatp2;
matMaster(:,:,3)=xMatp;

% Determine maximums and minimums
maxNoise=max(max(xMatNoise));
maxp2=max(max(xMatp2));
maxp=max(max(xMatp));

minNoise=min(min(xMatNoise));
minp2=min(min(xMatp2));
minp=min(min(xMatp));

MAX=max([maxNoise,maxp,maxp2]);
MIN=min([minNoise,minp,minp2]);

% Initialize standard deviation vector. NOTE: This vector is reset at each
% timestep.
SDVec=zeros(1,totCells-2);

% Hardcoded 
xyaxis=[1-10 totCells+10 -450 450];

figure

for i=start:fin
    
    % True solution
    subplot(2,2,1)
    plot(matMaster(i,:,1));
    title(['True: Step = ' num2str(i) ', Time (s) = '...
        num2str(measTime(i))]);
    xlabel('Cell');
    ylabel('Density');
    axis(xyaxis);
    
    % Est solution
    subplot(2,2,2)
    plot(matMaster(i,:,2));
    title(['Open Loop (No Data): Step = ' num2str(i) ', Time (s) = '...
        num2str(measTime(i))]);
    xlabel('Cell');
    ylabel('Density');
    axis(xyaxis);
    
    % Data solution-posterior
    subplot(2,2,3)
    plot(matMaster(i,:,3));
    hold on
    
    % Plot covariance bands with +/- std dev on posterior plot
    
    for j=1:totCells-2
        SDVec(1,j)=sqrt(Pk_p(j,j,i));
    end
    
    plot(matMaster(i,:,3)+sd*SDVec,'r');
    plot(matMaster(i,:,3)-sd*SDVec,'r');
    
    title(['Data (Posterior): Step = ' num2str(i) ', Time (s) = '...
        num2str(measTime(i))]);
    xlabel('Cell');
    ylabel('Density');
    axis(xyaxis);
   
    hold off
    
    legend('Posterior',['-',num2str(sd),'\sigma'],['+',num2str(sd),'\sigma'],'Location','SouthWest');
    
    pause(pDuration);
    
end