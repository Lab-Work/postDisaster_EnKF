% plotResults_indiv: plots all results in different figures. 
%
% NOTE: Not a function. Requires main code to be run first. Not the actual
% plotting function, but calls on the function (densPlotter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% WITH EQ

% True solution
figure
densPlotter(xMatNoise,mapLinksNoise,numLinks,measTime);

% Open loop solution
figure
densPlotter(xMatp2_wEQ,mapLinksApp2_wEQ,numLinks,measTime);

% Filtered solution
figure
densPlotter(xMatm_wEQ,mapLinksApp_wEQ,numLinks,measTime);
figure
densPlotter(xMatp_wEQ,mapLinksApp_wEQ,numLinks,measTime);

% WITH NO EQ

% True solution
figure
densPlotter(xMatNoise,mapLinksNoise,numLinks,measTime);

% Open loop solution
figure
densPlotter(xMatp2_nEQ,mapLinksApp2_nEQ,numLinks,measTime);

% Filtered solution
figure
densPlotter(xMatm_nEQ,mapLinksApp_nEQ,numLinks,measTime);
figure
densPlotter(xMatp_nEQ,mapLinksApp_nEQ,numLinks,measTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extremum check

% Add (:,1:end-2) to reshape the matrices in such a way to omit the boundary cells
extremumCheck(xMatNoise,xMatp2_wEQ,xMatm_wEQ,xMatp_wEQ)
extremumCheck(xMatNoise,xMatp2_nEQ,xMatm_nEQ,xMatp_nEQ)