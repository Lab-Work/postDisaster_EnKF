% setInitialParams: set the jam desnity and max flow to the mapLink objects
% before the filter runs
%
% INPUTS:
% mapLinks: a mapLinks object
% numLinks: integer for the number of links in the system

function mapLinks=setInitialParams(mapLinks,numLinks)

keyLinks=keys(mapLinks);

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    % Get necessary params
    vmax=clink.vmax;
    
    % Get the current capacity
    clink.qmax=clink.maxCap;
    qmax=clink.qmax;
    
    clink.rhoj=clink.maxRhoj;
    rhoj=clink.rhoj;
    
    % Get w
    clink.w=qmax/(rhoj-qmax/vmax);
    
    % Overwrite link object
    mapLinks(ckey)=clink;
    
end
    
    