% plotQmax: plots the distribution of qmax at the desired time steps. the
% left plot shows the true qmax values for each link. The middle plot shows
% the open loop solution distribution and the right plot shows the filter
% solution distribution
%
% INPUTS
% keyLinks: a keyLinks obj
% mapLinksNoise: the true mapLinks obj
% qmaxArray: array containing qmax values for all ensembles at all times
% (filter)
% qmaxArray2: array containing qmax values for all ensembles at all times
% (open loop)
% mapLinksApp_array: array containing mapLinks for all ensembles at all
% times (filter)
% mapLinksApp2_array: array containing mapLinks for all ensembles at all
% times (open loop)
% start: start time step
% fin: end time step
% noEn: integer for the number of ensembles
% pDuration: pause duration (may or may not be used) 
%
% NOTE: several parts are currently hardcoded 

function plotQmax(keyLinks,mapLinksNoise,qmaxArray,qmaxArray2,...
    mapLinksApp_array,mapLinksApp2_array,start,fin,noEn,...
    measTime,pDuration)

% FOR TRUE MODEL

numLinks=length(mapLinksNoise);

% Initialize
linkVec=zeros(1,numLinks);
qmaxVecAcc=zeros(1,numLinks);

% Create title
graphTitle='Link Order: ';

for i=1:numLinks
    
    % Create x axis
    linkVec(i)=i;
    
    % Current key
    ckey=keyLinks{i};
    
    % Current link
    clink=mapLinksNoise(ckey);
    
    % Extract qmax
    qmaxVecAcc(i)=clink.qmax;
    
    % Make title
    if i==numLinks
        graphTitle=[graphTitle,clink.ID];
    else
        graphTitle=[graphTitle,clink.ID,', '];
    end
    
end

figure 

% Plot
subplot(1,3,1)
plot(linkVec,qmaxVecAcc,'*');
title(graphTitle)
axis([0 numLinks+1 0 2500]); % hardcoded
xlabel('Link')
ylabel('qmax')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Qmax histograms for bridges

% A DAMAGE MAT NEEDS TO BE CREATED FOR EACH OBJ. THIS IS HARDCODED FOR NOW
% FOR THE SINGLE BRIDGE SYSTEM
damageMat=zeros(4,length(start:fin)); 
linkInd=2; % hardcoded for link three
maxCap=2000; % hardcoded

% Create a matrix that knows the value of each ensemble 
sampleMat=zeros(noEn,length(start:fin));

for k=start:fin
    
    % Reset counts
    INcount=0;
    MOcount=0;
    HVcount=0;
    CPcount=0;
    
    for i=1:noEn
        ckey=keyLinks{linkInd};
        maxCap=mapLinksApp_array{i}(ckey).maxCap;
        
        quant=qmaxArray(linkInd,1,i,k);
        
        if quant==maxCap
            INcount=INcount+1;
        elseif quant==maxCap*.75
            MOcount=MOcount+1;
        elseif quant==maxCap*.5
            HVcount=HVcount+1;
        else
            CPcount=CPcount+1;
        end
        
        % Keep track of all qmax values
        sampleMat(i,k)=quant;
        
    end
    
    % Keep track of amount of occurences for histogram
    damageMat(1,k)=INcount;
    damageMat(2,k)=MOcount;
    damageMat(3,k)=HVcount;
    damageMat(4,k)=CPcount;
    
end

% A DAMAGE MAT NEEDS TO BE CREATED FOR EACH OBJ
damageMat2=zeros(4,length(start:fin)); 

% Create a matrix that knows the value of each ensemble 
sampleMat2=zeros(noEn,length(start:fin));

for k=start:fin
    
    % Reset counts
    INcount=0;
    MOcount=0;
    HVcount=0;
    CPcount=0;
    
    for i=1:noEn
        ckey=keyLinks{linkInd};
        maxCap=mapLinksApp2_array{i}(ckey).maxCap;
        
        quant=qmaxArray2(linkInd,1,i,k);
        
        if quant==maxCap
            INcount=INcount+1;
        elseif quant==maxCap*.75
            MOcount=MOcount+1;
        elseif quant==maxCap*.5
            HVcount=HVcount+1;
        else
            CPcount=CPcount+1;
        end
        
        % Keep track of all qmax values
        sampleMat2(i,k)=quant;
        
    end
    
    % Keep track of amount of occurences for histogram
    damageMat2(1,k)=INcount;
    damageMat2(2,k)=MOcount;
    damageMat2(3,k)=HVcount;
    damageMat2(4,k)=CPcount;
    
end

% Compute mean and variance of qmax at each time step
meanVec=zeros(length(start:fin));
varVec=zeros(length(start:fin));

meanVec2=zeros(length(start:fin));
varVec2=zeros(length(start:fin));

count=1;

for k=start:fin
    
    % Use Matlabs built in functions
    meanVec(count)=mean(sampleMat(:,k));
    varVec(count)=var(sampleMat(:,k));
    
    meanVec2(count)=mean(sampleMat2(:,k));
    varVec2(count)=var(sampleMat2(:,k));
    
    count=count+1;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots

count=1;

for k=start:fin
   
    % Open loop 
    
    subplot(1,3,2)
    bar(1:4,damageMat2(:,k));
    title(['Forward, t=',num2str(measTime(k)),', \mu=',num2str(meanVec2(count)),...
        ', \sigma=',num2str(sqrt(varVec2(count)))])
    
    set(gca,'YLim',[0 100])
    set(gca,'XTickLabel',{'IN' 'MD' 'HI' 'CP'})
    
    xlabel('Damage state dist')
    ylabel(['Number of ens at t=',num2str(measTime(k))])
    
    % Kalman filter
    
    subplot(1,3,3)
    bar(1:4,damageMat(:,k));
    title(['Kalman, t=',num2str(measTime(k)),', \mu=',num2str(meanVec(count)),...
        ', \sigma=',num2str(sqrt(varVec(count)))])
    
    set(gca,'YLim',[0 100])
    set(gca,'XTickLabel',{'IN' 'MD' 'HI' 'CP'})
    
    xlabel('Damage state dist')
    ylabel(['Number of ens at t=',num2str(measTime(k))])
    
    pause(pDuration)
 
        % Manually operate pause
%     pause
    
    count=count+1;
    
end