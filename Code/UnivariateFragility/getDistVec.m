% getDistVec: takes in the mapLinks object and gives a column vector of the
% distances of the bridges from the EQ. The number of entries is equivalent
% to the number of bridges in the system.
%
% INPUTS
% mapLinks: a mapLinks object
% keyLinks: a cell array of link IDs
% numBridges: integer for the number of bridges in the system
% numLinks: integer for the number of links in the system

function distVec=getDistVec(mapLinks,keyLinks,numBridges,numLinks)

% Initialize the vector
distVec=zeros(numBridges,1);

% Keep track of the distVec
ind=1;

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    % Change bridges only
    if clink.isBridge~=0
        
        % Populate and increment counter
        distVec(ind)=clink.distFromEQ;
        ind=ind+1;
        
    end

end



