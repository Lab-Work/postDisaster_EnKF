% addInitialCond: adds the initial condition into the map obj
%
% INPUTS
% mapLinks: a mapLinks object
% keyLinks: a cell array of link IDs
% x0: vector which has the initial conditions for the links

function mapLinks=addInitialCond(mapLinks,keyLinks,x0)

for i=1:length(x0)
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    % Add initial density
    clink.rho0=x0(i);
    
    % Set initial qmax to the maximum capacity
    clink.qmax=clink.maxCap;
    
    % Overwrite link object
    mapLinks(ckey)=clink;
    
end
    