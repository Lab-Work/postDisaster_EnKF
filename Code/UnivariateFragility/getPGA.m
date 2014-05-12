% getPGA: takes in magnitude, distance, and site parameters and gives a
% vector of PGA values as dervied using the Campbell 1997 attenuation law.
%
% INPUTS
% M: an integer for the magnitude
% R: a vector of distances
% F: an integer of the fault type
% keyLinks: a cell array of link IDs
% mapBridges: a mapBridges object

function PGA=getPGA(M,R,F,mapLinks,keyLinks,mapBridges)

% Initialize
PGA=zeros(length(R),1);

numLinks=length(keyLinks);

% Keep track of the index
ind=1;

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    % Change bridges only
    if clink.isBridge~=0
        
        bkey=clink.isBridge;
        cbridge=mapBridges(bkey);
        
        % Get site condition parameters
        S_SR=cbridge.ssr;
        S_HR=cbridge.shr;
        
        % Campbell 1997 attentuation law
        lnPGA(ind)=-3.512+0.904*M-1.328*log(sqrt(R(ind)^2+(0.149*exp(0.647*M))^2))+...
            (1.125-0.112*log(R(ind))-0.0957*M)*F+(0.440-0.171*log(R(ind)))*S_SR+...
            (0.405-0.222*log(R(ind)))*S_HR;
        
        % Convert lnPGA to PGA (units of g)
        PGA(ind)=exp(lnPGA(ind));
        
        % Increment the index
        ind=ind+1;
        
    end
    
end
