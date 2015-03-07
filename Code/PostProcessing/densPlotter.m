% densPlotter: plots the densities in time and space in one graph
%
% INPUTS:
% xMat: matrix of solution densities (row: time, col: cell)
% mapLinks: a mapLinks object
% numLinks: integer for number of links in system
% measTime: vector of time steps (in seconds)

function densPlotter(xMat,mapLinks,numLinks,measTime)

numCells=size(xMat,2);

keyLinks=keys(mapLinks);

% Compute the max jam density considering all links

rhojvec=zeros(1,length(keyLinks));

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    rhoj=clink.maxRhoj;
    noLanes=clink.noLanes;
    
    rhojvec(i)=rhoj*noLanes;
    
end

% Uncomment if you want to see the maximum
maxrhoj=max(rhojvec);

% Create x and t vectors (x and y axes)
x_vec=1:numCells-2;
t_vec=measTime;

% Take xMat in such a way to eliminate the ghost cells
imagesc(x_vec,t_vec,xMat(:,x_vec));

% imagesc flips the Y axis so we must flip it back
colorbar
caxis([0,maxrhoj]);
set(gca,'Ydir','normal');
title('Traffic Density (veh/km)');
xlabel('Cell');
ylabel('Time (sec)');