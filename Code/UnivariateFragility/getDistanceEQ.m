% getDistanceEQ: takes the location of the earthquake and writes a new map
% object. The order here that the mapLinks are changed is not important.
%
% INPUTS
% mapLinks: a mapLinks object
% numLinks: integer for the number of links in the system
% mapBridges: a mapBridges object

function mapLinks=getDistanceEQ(mapLinks,numLinks,mapBridges)

keyLinks=keys(mapLinks);

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    % Change bridges only
    if clink.isBridge~=0
        
        % Get the bridge
        bkey=clink.isBridge;
        cbridge=mapBridges(bkey);
        
        % Distance
        clink.distFromEQ=cbridge.distToEQ;
        
    end
    
    % Overwrite link object
    mapLinks(ckey)=clink;
    
end
        
        