% getLengthCoord: take node coordinates into account and calculate length
% as well as the xycoord of the middle of the link
%
% INPUTS
% mapLinks: a mapLinks object
% mapNodes: a mapNodes object
% numLinks: integer for the number of links

function mapLinks=getLengthCoord(mapLinks,mapNodes,numLinks)

% Get keys
keyLinks=keys(mapLinks);

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    % Split node IDs
    str_nodes=clink.nodes;
    nodeIDs=regexp(str_nodes,'\s+','split');
    
    % Get xy coordinates of nodes
    x1=mapNodes(nodeIDs{1}).xycoord(1);
    y1=mapNodes(nodeIDs{1}).xycoord(2);
    x2=mapNodes(nodeIDs{2}).xycoord(1);
    y2=mapNodes(nodeIDs{2}).xycoord(2);
    
    % Compute xy coord of link
    x=(x1+x2)/2;
    y=(y1+y2)/2;
    clink.xycoord=[x y];
  
    % Compute length of link
    l=sqrt((x2-x1)^2+(y2-y1)^2);
    clink.roadLength=l;
    
    % Rewrite link object
    mapLinks(ckey)=clink;
    
end