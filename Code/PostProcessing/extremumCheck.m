% extremumCheck: Gives the min and max of the densities considering the
% time and space axes.
%
% xMatNoise: true solution 
% xMatp2: open loop solution
% xMatm: filtered solution (prior)
% xMatp: filtered solution (posterior)

function extremumCheck(xMatNoise,xMatp2,xMatm,xMatp)

% Cut out upstream and downstream ghost cells
min_xMatNoise=min(min(xMatNoise(:,1:end-2)));
max_xMatNoise=max(max(xMatNoise(:,1:end-2)));

min_xMatp2=min(min(xMatp2(:,1:end-2)));
max_xMatp2=max(max(xMatp2(:,1:end-2)));

min_xMatm=min(min(xMatm(:,1:end-2)));
max_xMatm=max(max(xMatm(:,1:end-2)));

min_xMatp=min(min(xMatp(:,1:end-2)));
max_xMatp=max(max(xMatp(:,1:end-2)));

disp(['xMatNoise: min=' num2str(min_xMatNoise) ', max=' num2str(max_xMatNoise)]); 
disp(['xMatp2: min=' num2str(min_xMatp2) ', max=' num2str(max_xMatp2)]);
disp(['xMatm: min=' num2str(min_xMatm) ', max=' num2str(max_xMatm)]); 
disp(['xMatp: min=' num2str(min_xMatp) ', max=' num2str(max_xMatp)]);