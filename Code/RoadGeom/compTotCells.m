% compTotCells: gives the total number of cells including the 2 ghost cells
%
% INPUTS
% mapLinks: a mapLinks object
% numLinks: integer for the number of links

function tot=compTotCells(mapLinks,numLinks)

keyLinks=keys(mapLinks);

% Initialize

tot=0;

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    % Add
    tot=tot+clink.noCells;
    
end