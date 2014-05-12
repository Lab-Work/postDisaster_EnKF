% plotError: plots the error (true/prior),(true/posterior) as a contour plot
%
% INPUTS
% xMatTrue: matrix of true solution densities (row: time, col: cell)
% xMatp2: matrix of filtered solution (prior) densities (row: time, col: cell)
% xMatp: matrix of filtered solution (posterior) densities (row: time, col: cell)
% totCells: integer for the number of cells in system
% measTime: vector of time steps (in seconds)

function plotError(xMatTrue,xMatp2,xMatp,totCells,measTime)

% Redefine to get rid of params
xMatTrue=xMatTrue(:,1:totCells-2);
xMatp2=xMatp2(:,1:totCells-2);
xMatp=xMatp(:,1:totCells-2);

% Caclulate error matrices
errPriorTrue=xMatp2-xMatTrue;
maxErrPriorTrue=max(max(errPriorTrue));
minErrPriorTrue=min(min(errPriorTrue));
errPostTrue=xMatp-xMatTrue;
maxErrPostTrue=max(max(errPostTrue));
minErrPostTrue=min(min(errPostTrue));

% New figure
figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data/True

subplot(1,3,1)

% Create x and t vectors (x and y axes). The "-2" neglects the boundary
% cells.
x_vec=1:totCells-2;
t_vec=measTime;

% Take xMat in such a way to eliminate the ghost cells
imagesc(x_vec,t_vec,errPriorTrue(:,x_vec));

% imagesc flips the Y axis so we must flip it back
colorbar
set(gca,'Ydir','normal');
caxis([min(minErrPriorTrue,minErrPostTrue),max(maxErrPriorTrue,maxErrPostTrue)]);
title('Density Error (Prior-True)');
xlabel('Cell');
ylabel('Time [sec]');
set(gca,'fontsize',10)
set(findall(gcf,'type','text'),'fontsize',10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data/True

subplot(1,3,2)

% Take xMat in such a way to eliminate the ghost cells
imagesc(x_vec,t_vec,errPostTrue(:,x_vec));

% imagesc flips the Y axis so we must flip it back
colorbar
set(gca,'Ydir','normal');
caxis([min(minErrPriorTrue,minErrPostTrue),max(maxErrPriorTrue,maxErrPostTrue)]);
title('Density Error (Posterior-True)');
xlabel('Cell');
ylabel('Time [sec]');
set(gca,'fontsize',10)
set(findall(gcf,'type','text'),'fontsize',10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Difference error 

diffErr=errPostTrue-errPriorTrue;

subplot(1,3,3)

% Take xMat in such a way to eliminate the ghost cells
imagesc(x_vec,t_vec,diffErr(:,x_vec));

% imagesc flips the Y axis so we must flip it back
colorbar
set(gca,'Ydir','normal');
title('Difference Error (Posterior-Prior)');
xlabel('Cell');
ylabel('Time [sec]');
set(gca,'fontsize',10)
set(findall(gcf,'type','text'),'fontsize',10);