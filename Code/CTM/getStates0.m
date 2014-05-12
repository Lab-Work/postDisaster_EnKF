% getStates: creates the initial density vector by extracting the initial
% condition from the mapLinks object
%
% INPUTS
% mapLinks: a mapLinks object
% totCells: integer for the total number of cells in system
% numLinks: integer for the number of links in the system

function x=getStates0(mapLinks,totCells,numLinks)

keyLinks=keys(mapLinks);

% Initialize
x=zeros(totCells,1);

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    % cstate can be a vector (should be able to anyway)
    cstate=clink.rho0;
    
    % Give the correct initial conditions to the correct cells
    x(clink.startCell:clink.endCell)=cstate;
    
end
    
    

