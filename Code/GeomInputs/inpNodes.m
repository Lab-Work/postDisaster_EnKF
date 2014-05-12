% inpNodes: Create nodes by giving them an ID and a location in [x y]
% coordinates
%
% INPUT
% numNodes: integer that tells how many nodes are in the network

function [valueNodes,keyNodes]=inpNodes(numNodes)

% Initialize
valueNodes=cell(1,numNodes);
keyNodes=cell(1,numNodes);

prompt={'ID','XY Coordinates (km)'};

for i=1:numNodes
    
    inp=inputdlg(prompt,['Node ' num2str(i)],1,{['n',num2str(i)],'0 0'});
    
    n=node;
    
    n.ID=inp{1};
    n.xycoord=str2num(inp{2});
    
    keyNodes{i}=n.ID;
    valueNodes{i}=n;
    
end
