% getStates0: creates a vector of the initial condition based on the
% intital conditions assigned to the links
%
% INPUTS
% mapLinks: a mapLinks object
% totCells: integer for number of cells in system
% numLinks: integer for number links in the system

function x=getStates0(mapLinks,totCells,numLinks)

keyLinks=keys(mapLinks);

% Initialize
x=zeros(totCells,1);

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    % cstate can be a vector 
    cstate=clink.rho0;
    
    % Assign initial conditions to appropriate cells
    x(clink.startCell:clink.endCell)=cstate;
    
end
    
    

