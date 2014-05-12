% distDraw0: draw a vector from a normal distribution with some mean vector
% and covariance matrix. Also the redraw count is returned.
%
% INPUTS
% meanVec: the initial denstities
% covMat: the initial covariance matrix
% mapLinks: a mapLinks object
% numLinks: integer for the number of links in the system
% totCells: integer for the number of cells in the system
% count: current redraw count
% enNo: the current ensemble number

function [drawvec count]=distDraw0(meanVec,covMat,mapLinks,numLinks,...
    totCells,count,enNo)

% Local count (refreshes per ensemble)
lCount=0;

% Initialization
drawvec=zeros(totCells,1);

keyLinks=keys(mapLinks);

for i=1:numLinks
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinks(ckey);
    
    % Get needed quantities
    startC=clink.startCell;
    endC=clink.endCell;
    rhoj=clink.rhoj;
    
    for j=startC:endC
        
        % Local cell count
        cCount=0;

        draw=normrnd(meanVec(j),sqrt(covMat(j,j)));
        
        % Redraw ensembles that are out of range
        while (draw<0||draw>rhoj)==1
            cCount=cCount+1;
            lCount=lCount+1;
            count=count+1;
            % Uncomment if you want to know which ensemnbles are initially
            % out of range
%             disp(['Redrawing initial ens: no: ' num2str(enNo)...
%                 ', redoNo: ' num2str(lCount) ', cellNo: ' num2str(j)]);
            draw=normrnd(meanVec(j),sqrt(covMat(j,j)));
        end
        
        % Populate
        drawvec(j)=draw;
        
    end
    
end