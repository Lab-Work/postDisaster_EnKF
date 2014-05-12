% updateLinks: uses the flows computed in updateNodes to run the CTM
%
% INPUTS
% xn: vector of densities at the current time step
% mapLinks: a mapLinks object
% numLinks: integer for the number of links in system
% totCells: integer for the number of cells in the system
% DT: time step (in seconds)
% leftBC: upstream ghost cell density
% rightBC: downstream ghost cell density
% isApp: integer that says whether or not this is an approximation (open
% loop or filter) (N=0,Y=1)
% noiseBC: vector containing the uncertainty on the boundary conditions
%
% The following are inputs but not used in the deterministic model
% 
% err_S: mean error in sending function
% err_Q: mean error in max flow region
% err_R: mean error in receiving function
% Q_S: variance in sending function
% Q_Q: variance in max flow region
% Q_R: variance in receiving function

function xn1=updateLinks(xn,mapLinks,numLinks,totCells,DT,leftBC,rightBC,...
    isApp,err_S,err_Q,err_R,Q_S,Q_Q,Q_R,noiseBC)

xn1=zeros(1,totCells);

keyLinks=keys(mapLinks);

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    % Get desired indices
    inds=clink.startCell:clink.endCell;
        
    % Get desired portion of vector to propagate through CTM
    xnlink=xn(inds);
    
    % Internal links
    if strcmp(clink.type,'i')==1
        
        % Run CTM
        xn1link=CTMnet(xnlink,clink,DT,isApp,err_S,err_Q,err_R,Q_S,Q_Q,Q_R);
        
    % Boundary (ghost cell) links
    else
        
        % Upstream boundary
        if strcmp(clink.type,'s')==1
            
            if isApp==1
                value=leftBC+normrnd(0,(noiseBC(1)/100)*leftBC);
            else
                value=leftBC;
            end
        
        % Downstream boundary    
        else
            
            if isApp==1
                value=rightBC+normrnd(0,(noiseBC(2)/100)*rightBC);
            else
                value=rightBC;
            end
            
        end
        
        xn1link=updateBoundary(value);
        
    end
    
    xn1(inds)=xn1link;
end

    