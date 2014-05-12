% updateNodes: updates the trueInflow and trueOutflow parameters of the
% links objects.
%
% INPUTS
% mapNodes: a mapNodes object
% numNodes: integer for the number of nodes in system
% mapLinks: a mapLinks object
% numLinks: integer for the number of links in system
% x: vector of densities at teh current time step
% isApp: integer that says whether or not this is an approximation (open
% loop or filter) (N=0,Y=1)
%
% The following are inputs but not used in the deterministic model
% 
% err_S: mean error in sending function
% err_Q: mean error in max flow region
% err_R: mean error in receiving function
% Q_S: variance in sending function
% Q_Q: variance in max flow region
% Q_R: variance in receiving function
%
% NOTE: this is able to handle merges and diverges (2 links -->1 link or
% 1 link -->2 links)

function mapLinks=updateNodes(mapNodes,numNodes,mapLinks,numLinks,x,isApp,...
    err_S,err_Q,err_R,Q_S,Q_Q,Q_R)

keyNodes=keys(mapNodes);
keyLinks=keys(mapLinks);

% For parameter estimation using a triangular model, need to express w in terms
% of qmax, vmax, and rhoj before any updating of nodes occurs

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    % Get necessary link params
    vmax=clink.vmax;
    qmax=clink.qmax;
    rhoj=clink.rhoj;
    
    % Compute w based on triangular model
    clink.w=qmax/(rhoj-qmax/vmax);
    
    % Overwrite link object
    mapLinks(ckey)=clink;
    
end

for i=1:numNodes
    
    % Current key
    ckey=keyNodes{i};
    
    % Current node
    cnode=mapNodes(ckey);
    
    % Split the links in and links out 
    IDs_in=regexp(cnode.links_in,'\W+','split');
    IDs_out=regexp(cnode.links_out,'\W+','split');
    
    % Merge
    if (length(IDs_in)>length(IDs_out))
        
        [nLinks,lkeys]=updateMergeQ(mapLinks,cnode,x,isApp,...
            err_S,err_Q,err_R,Q_S,Q_Q,Q_R);
        
    % Diverge
    elseif (length(IDs_in)<length(IDs_out))
        
        [nLinks,lkeys]=updateDivergeQ(mapLinks,cnode,x,isApp,...
            err_S,err_Q,err_R,Q_S,Q_Q,Q_R);
        
    % Normal Transition    
    else
        
        [nLinks,lkeys]=updateQ(mapLinks,cnode,x,isApp,...
            err_S,err_Q,err_R,Q_S,Q_Q,Q_R);
       
    end
    
    % Overwrite previous link objects in map with new ones
   
    for j=1:length(lkeys)
        
        mapLinks(lkeys{j})=nLinks(j);
    
    end
    
end