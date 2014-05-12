% getBridges: determines the links that are bridges and overwrites map
% object
%
% INPUTS
% mapLinks: a mapLinks object
% numLinks: integer for the number of links
% mapBridges: a mapBridges object
% numBridges: integer for the number of bridges

function mapLinks=getBridges(mapLinks,numLinks,mapBridges,numBridges)

% Get keys
keyLinks=keys(mapLinks);
keyBridges=keys(mapBridges);

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    for j=1:numBridges
        
        % Get the bridge
        bkey=keyBridges{j};
        cbridge=mapBridges(bkey);
        
        if isequal(clink.ID,cbridge.attachedToLink)
            
            clink.isBridge=cbridge.ID;
            
        end
        
    end
    
    % Rewrite link object
    mapLinks(ckey)=clink;
    
end