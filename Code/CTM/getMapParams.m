% getMapParams: overwrites the old map object with the damaged (or
% undamaged) mapLinks object.
%
% INPUTS
% mapLinks: a mapLinks object
% numLinks: integer for the number of links in the system
% keyLinks: a cell array of link IDs
% DSlims: vector designating ranges of different damage states
%
% NOTE: This is different from getMapParamsN, because there is no longer a
% damaged state forced (i.e. this is a stochastic process)

function mapLinks=getMapParams(mapLinks,numLinks,keyLinks,DSlims)

% Keep track of the number of bridges
ind=1;

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    % Get the max capacity
    maxCap=clink.maxCap;
    maxRhoj=clink.maxRhoj;
    
    % Check to see if the current ink is a bridge
    if clink.isBridge~=0
        
        % Uniform distribution from 0 to 1
        draw=rand;
        
        % Determine the damage state for that draw
        [ratio,DState]=compCap(DSlims(ind,:),draw);
        
        % Assign to clink
        clink.qmax=ratio*maxCap;
        clink.rhoj=ratio*maxRhoj;
        
        % Increment the index
        ind=ind+1;
        
    else % not a bridge, and thus no reduction in capacity
        
        clink.qmax=maxCap;
        clink.rhoj=maxRhoj;
        
    end
    
    % Overwrite old link object
    mapLinks(ckey)=clink;
    
end
