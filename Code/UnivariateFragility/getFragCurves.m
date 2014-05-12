% getFragCurves: takes in the median PGA and dispersion parameters for
% different limit states and creates fragility curves
%
% INPUTS
%
% PGAVec: vector of PGA values for fragility model (defines range)
% mapLinks: a mapLinks object
% keyLinks: a cell array of link IDs
% mapBridges: a mapBridges object
% numBridges: integer for the number of bridges in the system

function p=getFragCurves(PGAVec,mapLinks,keyLinks,mapBridges,numBridges)

% Initialize
numCurves=3; % hardcoded
p=zeros(numCurves,length(PGAVec),numBridges);

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
        
        % Get the bridge
        bkey=clink.isBridge;
        cbridge=mapBridges(bkey);
        
        % Get parameters
        med=cbridge.med;
        dsp=cbridge.dsp;
            
        for j=1:numCurves
            
            p(j,:,ind)=normcdf(log(PGAVec),log(med(j)),dsp);
            
        end
        
        % Increment
        ind=ind+1;
        
    end
    
end
        

%%%%%%%%%%%%%        
        
% % Initialize
% PGA=zeros(length(R),1);
% 
% numLinks=length(keyLinks);
% 
% % Keep track of the index
% ind=1;
% 
% for i=1:numLinks
%     
%     % Current key
%     ckey=keyLinks{i};
%     
%     % Current link
%     clink=mapLinks(ckey);
%     
%     % Change bridges only
%     if clink.isBridge~=0
%         
%         bkey=clink.isBridge;
%         cbridge=mapBridges(bkey);
%         
%         S_SR=cbridge.ssr;
%         S_HR=cbridge.shr;
%         
%         lnPGA(ind)=-3.512+0.904*M-1.328*log(sqrt(R(ind)^2+(0.149*exp(0.647*M))^2))+...
%             (1.125-0.112*log(R(ind))-0.0957)*F+(0.440-0.171*log(R(ind)))*S_SR+...
%             (0.405-0.222*log(R(ind)))*S_HR;
%         
%         % Convert lnPGA to PGA (units of g)
%         PGA(ind)=exp(lnPGA(ind));
%         
%         ind=ind+1;
%         
%     end
%     
% end

%%%%%%%

% p=zeros(length(med),length(PGAVec));
% 
% % Fragility curves 
% for i=1:length(med)
%         
%         p(i,:)=normcdf(log(PGAVec),log(med(i)),dsp(i));
%         
% end