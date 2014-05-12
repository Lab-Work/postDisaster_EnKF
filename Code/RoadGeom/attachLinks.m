% attachLinks: attach the incoming and outgoing links onto the correct
% nodes in map obj
% 
% INPUTS
% mapNodes: a mapNodes object
% numNodes: integer for the number of nodes
% mapLinks: a mapLinks object
% numLinks: integer for the number of links

function mapNodes=attachLinks(mapNodes,numNodes,mapLinks,numLinks)

keyNodes=keys(mapNodes);

for i=1:numNodes
    
    % Character counts
    inCount=1;
    outCount=1;
    
    % Current key
    ckeyN=keyNodes{i};
    
    % Current link
    cnode=mapNodes(ckeyN);
    
    keyLinks=keys(mapLinks);
    
    for j=1:numLinks
        
        % Current key
        ckeyL=keyLinks{j};
        
        % Current link
        clink=mapLinks(ckeyL);
        
        % Split node IDs
        str_nodes=clink.nodes;
        nodeIDs=regexp(str_nodes,'\s+','split');
                
        % Separate the boundaries. The awkward bookkeeping is necessary
        % for merges and diverges
        if strcmp(clink.type,'s')==1
            
            if strcmp(nodeIDs{2},cnode.ID)==1
                cnode.links_in(inCount:inCount+length(clink.ID)-1)=clink.ID;
                inCount=inCount+length(clink.ID)+1;
            end
            
        elseif strcmp(clink.type,'e')==1
            
            if strcmp(nodeIDs{1},cnode.ID)==1
                cnode.links_out(outCount:outCount+length(clink.ID)-1)=clink.ID;
                outCount=outCount+length(clink.ID)+1;
            end
          
        else
            
            if strcmp(nodeIDs{2},cnode.ID)==1
                cnode.links_in(inCount:inCount+length(clink.ID)-1)=clink.ID;
                inCount=inCount+length(clink.ID)+1;
            end
            
            if strcmp(nodeIDs{1},cnode.ID)==1
                cnode.links_out(outCount:outCount+length(clink.ID)-1)=clink.ID;
                outCount=outCount+length(clink.ID)+1;
            end
            
        end
        
        % Rewrite node object
        mapNodes(ckeyN)=cnode;
        
    end
    
end

