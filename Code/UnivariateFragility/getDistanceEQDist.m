% getDistanceEQDist: takes the distance of the earthquake and writes a new map
% object. The order here that the mapLinks are changed is not important.
% The distance is drawn from a normal distribution centered around the true
% distance
%
% INPUTS
% mapLinks: a mapLinks object
% numLinks: integer for the number of links in the system
% mapBridges: a mapBridges object
% isEQInp: integer that tells the model whether to consider the EQ or not

function mapLinks=getDistanceEQDist(mapLinks,numLinks,mapBridges,isEQInp)

keyLinks=keys(mapLinks);

% Percent standard deviation noise
noiseDist=20;

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
        
        % Draw the distance from a normal distribution
        if isEQInp==1
            clink.distFromEQ=normrnd(cbridge.distToEQ,(noiseDist/100)*...
                cbridge.distToEQ);
        else
            clink.distFromEQ=1e6;
        end
        
    end
    
    % Overwrite the link object
    mapLinks(ckey)=clink;
    
end
        
        