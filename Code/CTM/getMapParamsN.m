% getMapParamsN: overwrites the old mapLinks object with the damaged
% mapLinks object. THIS IS A SPECIAL INSTANCE OF THE FUNCTION WHERE WE
% WANT THE BRIDGE TO BE DAMAGED!
%
% INPUTS
% mapLinks: a mapLinks object
% numLinks: integer for the number of links in the system
% keyLinks: a cell array of link IDs
% mapBridges: a mapBridges object
%
% NOTE: In the true model, we force a damage state, depending on the type of
% earthquake object created (i.e. it is deterministic)

function mapLinks=getMapParamsN(mapLinks,numLinks,keyLinks,mapBridges)

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
    
    % Check to see if the current link is a bridge
    if clink.isBridge~=0
        
        bkey=clink.isBridge;
        cbridge=mapBridges(bkey);
        
        % Prescribe damage for true solution based on the bridge
        if isequal(cbridge.trueDmg,'MD')
            linkCap=0.75*maxCap;
            linkRhoj=0.75*maxRhoj;
        elseif isequal(cbridge.trueDmg,'HI')
            linkCap=0.5*maxCap;
            linkRhoj=0.5*maxRhoj;
            elseif isequal(cbridge.trueDmg,'TO')
            linkCap=0*maxCap;
            linkRhoj=0*maxRhoj;
        else 
            linkCap=maxCap;
            linkRhoj=maxRhoj;
        end
        
        % Display text
        disp(['Link: ' clink.ID ', Capacity: ' num2str(linkCap)]);
        
        % Increment the index
        ind=ind+1;
        
    else % not a bridge, and thus no reduction in capacity
        
        linkCap=maxCap;
        linkRhoj=maxRhoj;
        
    end
    
    clink.qmax=linkCap;
    clink.rhoj=linkRhoj;
    
    % Overwrite old link object
    mapLinks(ckey)=clink;
    
end
