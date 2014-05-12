% plotResults: plots all results in different figures. EQ or no EQ results
% only (not both)
%
% NOTE: Not a function. Requires main code to be run first. Not the actual
% plotting function, but calls on the function (densPlotter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% True solution
figure
densPlotter(xMatNoise,mapLinksNoise,numLinks,measTime);

% Open loop solution
figure
densPlotter(xMatp2,mapLinksApp2,numLinks,measTime);

% Filtered solution
figure
densPlotter(xMatm,mapLinksApp,numLinks,measTime);
figure
densPlotter(xMatp,mapLinksApp,numLinks,measTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extremum check

% Add (:,1:end-2) to reshape the matrices in such a way to omit the boundary cells
extremumCheck(xMatNoise,xMatp2,xMatm,xMatp)