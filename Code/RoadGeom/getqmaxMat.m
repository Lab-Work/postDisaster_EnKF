% getqmaxMat: create a vector of the qmax value for each link
%
% INPUTS
% mapLinks: a mapLinks object
% keyLinks: a cell array of link IDs

function qmaxMat=getqmaxMat(mapLinks,keyLinks)

numLinks=length(mapLinks);

% Initialize
qmaxMat=zeros(numLinks,1);

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    % Extract qmax
    qmaxMat(i)=clink.qmax;
    
end



