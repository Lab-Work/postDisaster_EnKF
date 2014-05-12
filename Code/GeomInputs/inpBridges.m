% inpBridges: Give the IDs of the links which you want to be bridges.
% Also describes the bridge material and site conditions.
% 
% INPUTS
% numBridges: integer that tells how many bridges are in the network

function [valueBridges,keyBridges]=inpBridges(numBridges)

% Initialize
valueBridges=cell(1,numBridges);
keyBridges=cell(1,numBridges);

% prompt needs ID,type,conn nodes,free flow speed,backward propagating wave
% speed, maximum flow capacity, number of lanes

prompt={'ID','Associated link ID','Ssr (soil)',...
    'Shr (soil)','Fragility (from Evensen)','Dist to EQ (km)',...
    'Dmg in True Mdl (IN,MD,HI,TO)'};

for i=1:numBridges
    
    inp=inputdlg(prompt,['Bridge ' num2str(i)],1,...
        {['b',num2str(i)],'link ID','0','0','MSC_S','0','damage'});
    
    b=bridge;
    
    b.ID=inp{1};
    b.attachedToLink=inp{2};
    b.ssr=str2num(inp{3});
    b.shr=str2num(inp{4});
    b.type=inp{5};
    % NOTE: this is coded for only four types of bridges right now with 3
    % limit states 
    if isequal(inp{5},'MSC_C')
        b.med=[.15 .52 1.03];
        b.dsp=.70;
    elseif isequal(inp{5},'MSC_S')
        b.med=[.18 .31 .5];
        b.dsp=.55;
    elseif isequal(inp{5},'MSSS_C')
        b.med=[.20 .57 1.17];
        b.dsp=.65;
    else isequal(inp{5},'MSSS_S')
        b.med=[.24 .44 .82];
        b.dsp=.5;
    end

    b.distToEQ=str2num(inp{6});
    b.trueDmg=inp{7};   
    
    keyBridges{i}=b.ID;
    valueBridges{i}=b;
    
end

%%%%%%%%%%
% prompt={'Link IDs'};
% 
% % Input bridges
% inp=inputdlg(prompt,'Which links are bridges?',1);
% 
% % Split link IDs
% str_bridges=inp{1};
% bridgeLinks=regexp(str_bridges,'\s+','split');


