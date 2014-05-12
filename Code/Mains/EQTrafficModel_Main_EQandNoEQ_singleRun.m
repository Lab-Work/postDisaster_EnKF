% EQTrafficModel_Main_EQandNoEQ_singleRun: the main file to 
% execute. This will run the same initial conditions with and without the 
% earthquake.
%
% In this script, variables with a "_wEQ" ending mean they are for the
% model that includes knowledge of the earthquake and variables with a
% "_nEQ" ending mean they are for the model that DOES NOT include knowledge
% of the earthquake.

tic

close all; clear all; clc

% Warnings:
warning('Triangular fundamental diagram is assumed');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % INPUT PARAMETERS (USER INPUTS)
% 
% % Nodes
% 
% % NOTE: mapNodes does not need to be copied so passing the map back through
% % the function should be fine, but pass back valueNodes just in case
% numNodes=4;
% [valueNodes,keyNodes]=inpNodes(numNodes);
% 
% % Links
% 
% % Pass back the cell of links instead of the actual map so multiple map
% % object with difrerent handles can be created
% numLinks=5; % internal/boundary combined total
% [valueLinks,keyLinks]=inpLinks(numLinks);
% 
% % Sensors
% 
% % % NOTE: mapSensors does not need to be copied so passing the map back through
% % % the function should be fine, but pass back valueSensors just in case
% % numSensors=18;
% % [valueSensors,keySensors]=inpSensors(numSensors);
% 
% % Tell which links are actually bridges
% numBridges=1;
% [valueBridges,keyBridges]=inpBridges(numBridges);
% 
% % EQ 
% newEQ=inpEQ;
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PRELOADED ROAD GEOMETRY

% Choose the desired damage state of the Caruthersville Bridge
dmgType=2; %1=MD, 2=HI, 3=TO

if dmgType==1
    load Geom_Midwest_SensNoise5_MD.mat
elseif dmgType==2
    load Geom_Midwest_SensNoise5_HI.mat
else
    load Geom_Midwest_SensNoise5_TO.mat
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SIMULATION AND FILTER VARIABLES

tmax=2400; % time of simulation [seconds]
DT=15; % time step [seconds]

% Obtain time vector
measTime=0:DT:tmax;

% Time at which EQ occurs (make sure it is a factor of tmax and a multiple of DT)
teq=600;
keq=(teq/DT)+1; % Kalman time step at which EQ occurs

% EnKF Variables
noEn=150; % number of ensemble members
tolerance=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BRIDGE VARIABLES (Nielson 2007) 

% Vector of PGA inputs to build fragility curve. These variables describe
% how accurately we want to map computed PGA values to our fragility curves
delta=.001; % The degree of accuracy
SF=1/delta; % Scale factor for rounding purposes
PGAVec=0:delta:3; % Choose a reasonable range

% Parameters for attenuation (Refer to Campbell 1997)
F=0; % Numerical value for type of fault, scalar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIAL AND BOUNDARY CONDITIONS, AND NOISE MODELS 

% Inital conditions. Currently the order is [internal links;left BC;right
% BC]

% Of accurate/noise solution (TOTAL NOT PER LANE)
rho0Acc=[30;30;30;40;20];

% Of approx (Kalman) solution (TOTAL NOT PER LANE)
rho0App=[10;10;10;30;10];

% Boundary conditions for acc and app models. This is dependent on the
% order of rho0App and rho0Acc
startAcc=rho0Acc(end-1); % left BC
endAcc=rho0Acc(end); % right BC
startApp=rho0App(end-1); % left BC
endApp=rho0App(end); % right BC

% Boundary condition noise (density percentages)
startSig=25;
endSig=25;
noiseBC=[startSig endSig];

% Normal noise parameters (mean, std dev)
err_S=0;
err_Q=0;
err_R=0;
Q_S=50^2;
Q_Q=100^2;
Q_R=150^2;

% Probability of sensor failure. Can also be set to a nonzero value, but
% this would change the H matrix.
pf_sens=0;

% Set whether or not the true model and approx model have imperfect params
% (1=Yes, 0=No)
isImperfParams=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CREATE MULTIPLE MAPLINKS/MAPNODES FOR EXACT/ACCURATE ETC

% For true solution (DT,rho0Acc)
mapLinksNoise=containers.Map(keyLinks,valueLinks);

% For Kalman solution (DT,rho0App). A second map object is necessary so that
% the open loop model can be run.

% Initialize
mapLinksApp_wEQ_array=cell(1,noEn);
mapLinksApp2_wEQ_array=cell(1,noEn);
mapLinksApp_nEQ_array=cell(1,noEn);
mapLinksApp2_nEQ_array=cell(1,noEn);

vmaxMat=zeros(noEn,numLinks);

for i=1:noEn
    
    if isImperfParams==1
        
        % For these solutions slightly perturb values of true parameters.
        [vmaxVec,nvalueLinks]=alterParams(valueLinks);
        mapLinksApp_wEQ_array{i}=containers.Map(keyLinks,nvalueLinks); % filter, EQ
        mapLinksApp2_wEQ_array{i}=containers.Map(keyLinks,nvalueLinks); % open loop, EQ
        mapLinksApp_nEQ_array{i}=containers.Map(keyLinks,nvalueLinks); % filter, no EQ
        mapLinksApp2_nEQ_array{i}=containers.Map(keyLinks,nvalueLinks); % open loop, no EQ
        
        vmaxMat(i,:)=vmaxVec;
        
    else
        
        % These are the same as the true model parameters
        mapLinksApp_wEQ_array{i}=containers.Map(keyLinks,valueLinks); % filter, EQ
        mapLinksApp2_wEQ_array{i}=containers.Map(keyLinks,valueLinks); % open loop, EQ
        mapLinksApp_nEQ_array{i}=containers.Map(keyLinks,valueLinks); % filter, no EQ
        mapLinksApp2_nEQ_array{i}=containers.Map(keyLinks,valueLinks); % open loop, no EQ
        
    end
    
end

% Get maximum of vmax for each link
maxvmaxVec=zeros(1,numLinks);

if isImperfParams==1
    for i=1:numLinks
        maxvmaxVec(i)=max(vmaxMat(:,i));
    end  
else
    for i=1:numLinks
        ckey=keyLinks{i};
        clink=mapLinksNoise(ckey);
        maxvmaxVec(i)=clink.vmax;
    end
end

% mapNodes
mapNodes=containers.Map(keyNodes,valueNodes);

% mapSensors (ONLY TRAFFIC DENSITY SENSORS)
mapSensors=containers.Map(keySensors,valueSensors);

% mapBridges
mapBridges=containers.Map(keyBridges,valueBridges);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE NECESSARY QUANTITIES FOR EACH MAPLINKS CONTAINER

% Compute length of each link and overwrite previous map obj.
mapLinksNoise=getLengthCoord(mapLinksNoise,mapNodes,numLinks);

for i=1:noEn
    
    % Unwrap
    mapLinksApp_wEQ=mapLinksApp_wEQ_array{i};
    mapLinksApp2_wEQ=mapLinksApp2_wEQ_array{i};
    mapLinksApp_nEQ=mapLinksApp_nEQ_array{i};
    mapLinksApp2_nEQ=mapLinksApp2_nEQ_array{i};
    
    mapLinksApp_wEQ=getLengthCoord(mapLinksApp_wEQ,mapNodes,numLinks);
    mapLinksApp2_wEQ=getLengthCoord(mapLinksApp2_wEQ,mapNodes,numLinks);
    mapLinksApp_nEQ=getLengthCoord(mapLinksApp_nEQ,mapNodes,numLinks);
    mapLinksApp2_nEQ=getLengthCoord(mapLinksApp2_nEQ,mapNodes,numLinks);
    
    % Rewrap
    mapLinksApp_wEQ_array{i}=mapLinksApp_wEQ;
    mapLinksApp2_wEQ_array{i}=mapLinksApp2_wEQ;
    mapLinksApp_nEQ_array{i}=mapLinksApp_nEQ;
    mapLinksApp2_nEQ_array{i}=mapLinksApp2_nEQ;
    
end

% Compute dx,number of cells, and start/end cell for link using CFL
% condition and overwrite previous map obj.
mapLinksNoise=getCellDiscImperf(DT,mapLinksNoise,keyLinks,numLinks,maxvmaxVec);

for i=1:noEn
    
    % Unwrap
    mapLinksApp_wEQ=mapLinksApp_wEQ_array{i};
    mapLinksApp2_wEQ=mapLinksApp2_wEQ_array{i};
    mapLinksApp_nEQ=mapLinksApp_nEQ_array{i};
    mapLinksApp2_nEQ=mapLinksApp2_nEQ_array{i};
    
    mapLinksApp_wEQ=getCellDiscImperf(DT,mapLinksApp_wEQ,keyLinks,numLinks,maxvmaxVec);
    mapLinksApp2_wEQ=getCellDiscImperf(DT,mapLinksApp2_wEQ,keyLinks,numLinks,maxvmaxVec);
    mapLinksApp_nEQ=getCellDiscImperf(DT,mapLinksApp_nEQ,keyLinks,numLinks,maxvmaxVec);
    mapLinksApp2_nEQ=getCellDiscImperf(DT,mapLinksApp2_nEQ,keyLinks,numLinks,maxvmaxVec);
    
    % Rewrap
    mapLinksApp_wEQ_array{i}=mapLinksApp_wEQ;
    mapLinksApp2_wEQ_array{i}=mapLinksApp2_wEQ;
    mapLinksApp_nEQ_array{i}=mapLinksApp_nEQ;
    mapLinksApp2_nEQ_array{i}=mapLinksApp2_nEQ;
    
end

% The total number of cells takes into account the CTM model and the two
% ghost cells.
totCells=compTotCells(mapLinksNoise,numLinks);

% Add the initial conditions and overwrite previous map obj
mapLinksNoise=addInitialCond(mapLinksNoise,keyLinks,rho0Acc);

for i=1:noEn
    
    % Unwrap
    mapLinksApp_wEQ=mapLinksApp_wEQ_array{i};
    mapLinksApp2_wEQ=mapLinksApp2_wEQ_array{i};
    mapLinksApp_nEQ=mapLinksApp_nEQ_array{i};
    mapLinksApp2_nEQ=mapLinksApp2_nEQ_array{i};
    
    mapLinksApp_wEQ=addInitialCond(mapLinksApp_wEQ,keyLinks,rho0App);
    mapLinksApp2_wEQ=addInitialCond(mapLinksApp2_wEQ,keyLinks,rho0App);
    mapLinksApp_nEQ=addInitialCond(mapLinksApp_nEQ,keyLinks,rho0App);
    mapLinksApp2_nEQ=addInitialCond(mapLinksApp2_nEQ,keyLinks,rho0App);
    
    % Rewrap
    mapLinksApp_wEQ_array{i}=mapLinksApp_wEQ;
    mapLinksApp2_wEQ_array{i}=mapLinksApp2_wEQ;
    mapLinksApp_nEQ_array{i}=mapLinksApp_nEQ;
    mapLinksApp2_nEQ_array{i}=mapLinksApp2_nEQ;
    
end

% Determine which links are bridges and overwrite previous map obj
mapLinksNoise=getBridges(mapLinksNoise,numLinks,mapBridges,numBridges);

for i=1:noEn
    
    % Unwrap
    mapLinksApp_wEQ=mapLinksApp_wEQ_array{i};
    mapLinksApp2_wEQ=mapLinksApp2_wEQ_array{i};
    mapLinksApp_nEQ=mapLinksApp_nEQ_array{i};
    mapLinksApp2_nEQ=mapLinksApp2_nEQ_array{i};
    
    mapLinksApp_wEQ=getBridges(mapLinksApp_wEQ,numLinks,mapBridges,numBridges);
    mapLinksApp2_wEQ=getBridges(mapLinksApp2_wEQ,numLinks,mapBridges,numBridges);
    mapLinksApp_nEQ=getBridges(mapLinksApp_nEQ,numLinks,mapBridges,numBridges);
    mapLinksApp2_nEQ=getBridges(mapLinksApp2_nEQ,numLinks,mapBridges,numBridges);
    
    % Rewrap
    mapLinksApp_wEQ_array{i}=mapLinksApp_wEQ;
    mapLinksApp2_wEQ_array{i}=mapLinksApp2_wEQ;
    mapLinksApp_nEQ_array{i}=mapLinksApp_nEQ;
    mapLinksApp2_nEQ_array{i}=mapLinksApp2_nEQ;
    
end

% Attach links onto nodes and overwrite previous map obj
mapNodes=attachLinks(mapNodes,numNodes,mapLinksApp_wEQ,numLinks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MEASUREMENT MATRIX, H

% Rearrange sensor order and give a global cell position. keySensors here
% gives the keys back in the desired order.
[mapSensors keySensors sortInd]=getSensorPos(mapSensors,mapLinksApp_wEQ,...
    numSensors);

% Attach sensors onto respective links. No need to attach sensors to open
% loop map objects.

for i=1:noEn
    
    % Unwrap
    mapLinksApp_wEQ=mapLinksApp_wEQ_array{i};
    mapLinksApp_nEQ=mapLinksApp_nEQ_array{i};
    
    mapLinksApp_wEQ=attachSensors(mapLinksApp_wEQ,mapSensors,numSensors,keySensors);
    mapLinksApp_nEQ=attachSensors(mapLinksApp_nEQ,mapSensors,numSensors,keySensors);
    
    % Rewrap
    mapLinksApp_wEQ_array{i}=mapLinksApp_wEQ;
    mapLinksApp_nEQ_array{i}=mapLinksApp_nEQ;
    
end

% Get measurement matrix, H and its dimensions, kk and n
[H kk n]=getHSens(mapSensors,keySensors,totCells,numSensors);
Htrue=H;

% Heq is the modified H matrix when the earthquake randomly damages some sensors
[Heq sensInd]=getHeq(mapSensors,keySensors,totCells,numSensors,pf_sens);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RUN CTM (forward progression) for true solution

% Initialize. NOTE: Parameters of the true model are given a "Noise" or "N"
% ending in their name.
xMatNoise=zeros(length(measTime),totCells);

% Initial conditions. Implementation depends on the mapSensors being crated
% in order of cell discretization
x0Noise=getStates0(mapLinksNoise,totCells,numLinks);

% For true model
xMatNoise(1,:)=x0Noise';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% THIS SECTION IS NOT NECESSARY, BUT COULD BE USEFUL FOR DIAGNOSTICS.

% For true model, determine qmax due to damage from the earthquake. 'N'
% here indicates those variables related with the true model.

% Compute bridge location w.r.t. epicenter
mapLinksNoise=getDistanceEQ(mapLinksNoise,numLinks,mapBridges);

% Construct vector of distances
distVecN=getDistVec(mapLinksNoise,keyLinks,numBridges,numLinks);

% Get peak ground accelerations
PGAN=getPGA(newEQ.mag,distVecN,F,mapLinksNoise,keyLinks,mapBridges);

% Create deterministic fragility curves for fragility model. 
p=getFragCurves(PGAVec,mapLinksNoise,keyLinks,mapBridges,numBridges);

% Takes the Sa from the attenuation law and gives the corresponding
% probability of exceedance
PexN=getPex(PGAN,p,SF,PGAVec);

% Take the probability of exceedance and compute the probability of damage
% states as a vector
ProbDSN=getProbVec(PexN);

% Get the range of the four limit states
DSlimsN=getDSRange(ProbDSN);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Forward CTM model
for k=1:length(measTime)-1
    
    % Assign the parameter based on the probability table. We only want to
    % assign ONCE at the instance the EQ occurs
    if k==keq
        
        mapLinksNoise=getMapParamsN(mapLinksNoise,numLinks,keyLinks,mapBridges);
        
    end
    
    % Density propogation via CTM
    xMatNoise(k+1,1:totCells)=forwardPredict(DT,xMatNoise(k,:),mapLinksNoise,...
        mapNodes,totCells,startAcc,endAcc,0,err_S,err_Q,err_R,Q_S,Q_Q,Q_R,noiseBC);
    
end

% Uncomment if you want to just run deterministic model
% densPlotter(xMatNoise,mapLinksNoise,numLinks,measTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% KALMAN FILTER PREPARATION

% Create a measurement noise matrix.
noisyMat=getNoisyMat(mapSensors,numSensors,keySensors,measTime);

% Get the measurement vector by adding the noisyMat to the true solution 
y=getyMeas(H,xMatNoise,measTime,numSensors,noisyMat);

% Kalman time
K=length(measTime);

% Initialize 3D arrays and 4D arrays that contain density values and
% covariance errors for all ensembles and all time steps (m=prior,
% p=posterior)
[xk_m_wEQ xk_p_wEQ Pk_m_wEQ Pk_p_wEQ]=initialize3DArray(n,K);
[xi_m_wEQ xi_p_wEQ]=initialize4DArray(n,noEn,K);

[xk_m_nEQ xk_p_nEQ Pk_m_nEQ Pk_p_nEQ]=initialize3DArray(n,K);
[xi_m_nEQ xi_p_nEQ]=initialize4DArray(n,noEn,K);

% qmaxArray to keep track of qmax of every ensemble for every link
qmaxArray_wEQ=zeros(numLinks,1,noEn,K);

qmaxArray_nEQ=zeros(numLinks,1,noEn,K);

% Averaged measurement noise
v_bar=zeros(numSensors,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GET INITIAL STATES AND INITIAL ERROR COV OF STATES. ONLY DO THIS PROCESS
% ONCE!

% Initial states
xk_p_wEQ(1:totCells,1,1)=getStates0(mapLinksApp_wEQ,totCells,numLinks);
xk_p_nEQ(1:totCells,1,1)=getStates0(mapLinksApp_nEQ,totCells,numLinks);

x0App_wEQ=xk_p_wEQ(:,1,1);
x0App_nEQ=xk_p_nEQ(:,1,1);

% Initial error covariance matrix
Pk_p_wEQ(:,:,1)=initialCov(x0Noise,x0App_wEQ,tolerance,totCells);
Pk_p_nEQ(:,:,1)=initialCov(x0Noise,x0App_nEQ,tolerance,totCells);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% THE SUBSCRIPT 2 IS FOR THE FORWARD PROPOGATION OF ENSEMBLES WITH NO
% SENSOR DATA (OPEN LOOP)

% Initialize (no need to create add'l model noise or meas noise arrays)
[xk_m2_wEQ xk_p2_wEQ Pk_m2_wEQ Pk_p2_wEQ]=initialize3DArray(n,K);

xk_p2_wEQ(:,1,1)=xk_p_wEQ(:,1,1);

[xi_m2_wEQ xi_p2_wEQ]=initialize4DArray(n,noEn,K);

qmaxArray2_wEQ=zeros(numLinks,1,noEn,K);

%%%

[xk_m2_nEQ xk_p2_nEQ Pk_m2_nEQ Pk_p2_nEQ]=initialize3DArray(n,K);

xk_p2_nEQ(:,1,1)=xk_p_nEQ(:,1,1);

[xi_m2_nEQ xi_p2_nEQ]=initialize4DArray(n,noEn,K);

qmaxArray2_nEQ=zeros(numLinks,1,noEn,K);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If you dont do this overwriting, qmax and rhoj will both be 0, which will
% cause distDraw0 to go in an infinite loop.
for i=1:noEn
    
    % Unwrap
    mapLinksApp_wEQ=mapLinksApp_wEQ_array{i};
    mapLinksApp2_wEQ=mapLinksApp2_wEQ_array{i};
    mapLinksApp_nEQ=mapLinksApp_nEQ_array{i};
    mapLinksApp2_nEQ=mapLinksApp2_nEQ_array{i};
    
    mapLinksApp_wEQ=setInitialParams(mapLinksApp_wEQ,numLinks);
    mapLinksApp2_wEQ=setInitialParams(mapLinksApp2_wEQ,numLinks);
    mapLinksApp_nEQ=setInitialParams(mapLinksApp_nEQ,numLinks);
    mapLinksApp2_nEQ=setInitialParams(mapLinksApp2_nEQ,numLinks);
    
    % Rewrap
    mapLinksApp_wEQ_array{i}=mapLinksApp_wEQ;
    mapLinksApp2_wEQ_array{i}=mapLinksApp2_wEQ;
    mapLinksApp_nEQ_array{i}=mapLinksApp_nEQ;
    mapLinksApp2_nEQ_array{i}=mapLinksApp2_nEQ;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RUN FILTER

redrawCount=0;

% Initial ensembles draws

for i=1:noEn
    
    % distDraw0 gives only the density draws
    [xi_p_wEQ(1:totCells,1,i,1) redrawCount]=distDraw0(xk_p_wEQ(:,1,1),Pk_p_wEQ(:,:,1),mapLinksApp_wEQ,...
        numLinks,totCells,redrawCount,i);
    
    % Open loop solution has same inital ensembles
    xi_p2_wEQ(:,1,i,1)=xi_p_wEQ(:,1,i,1);
    
    % No EQ solution has same initial ensembles
    xi_p_nEQ(:,1,i,1)=xi_p_wEQ(:,1,i,1);
    xi_p2_nEQ(:,1,i,1)=xi_p_wEQ(:,1,i,1);
    
end

redrawCount

% Initialize X5 and filtered ensembles arrays. This allows every ensemble
% to be reproduced at every step, even after the simulation has been
% performed.
X5_array_wEQ=cell(K,1);
Af_array_wEQ=cell(K,1);
D_array_wEQ=cell(K,1);
Gam_array_wEQ=cell(K,1);
innov_array_wEQ=cell(K,1);
kalGain_array_wEQ=cell(K,1);

X5_array_nEQ=cell(K,1);
Af_array_nEQ=cell(K,1);
D_array_nEQ=cell(K,1);
Gam_array_nEQ=cell(K,1);
innov_array_nEQ=cell(K,1);
kalGain_array_nEQ=cell(K,1);

% Initial filtered ensembles in array form. ONLY DO THIS ONCE!
Af_i=composeEnsToMat(xi_p_wEQ,totCells,noEn,1);
Af_array_wEQ{1}=Af_i;
Af_array_nEQ{1}=Af_i;

for k=1:K-1

    % If you want to make sure code is processing
    if mod(k,10)==0
        k
        toc
    end
    
    % MODEL PREDICTION (PRIOR)
    
    for i=1:noEn
        
        % Unwrap
        mapLinksApp_wEQ=mapLinksApp_wEQ_array{i};
        mapLinksApp2_wEQ=mapLinksApp2_wEQ_array{i};
        mapLinksApp_nEQ=mapLinksApp_nEQ_array{i};
        mapLinksApp2_nEQ=mapLinksApp2_nEQ_array{i};
        
        % FRAGILITY FUNCTIONS
        
        % If the EQ can occur
        if k>=keq
            
            % Get the EQ characteristics. Note that distribution parameters
            % are defined in the function itself. Only need one set of
            % params for wEQ and one set of params for nEQ. Open loop or
            % filter does not matter.
            
            m_wEQ=getEQParams(1,newEQ);
            m_nEQ=getEQParams(0,newEQ);
            
            % Create the earthquakes. Only need one set of
            % params for wEQ and one set of params for nEQ. Open loop or
            % filter does not matter.
            nEQ_wEQ=editEQ(m_wEQ);
            nEQ_nEQ=editEQ(m_nEQ);
            
            % Compute bridge location w.r.t. epicenter. Assign to all
            % mapLink objects.
            mapLinksApp_wEQ=getDistanceEQDist(mapLinksApp_wEQ,numLinks,mapBridges,1);
            mapLinksApp2_wEQ=getDistanceEQDist(mapLinksApp2_wEQ,numLinks,mapBridges,1);
            
            mapLinksApp_nEQ=getDistanceEQDist(mapLinksApp_nEQ,numLinks,mapBridges,0);
            mapLinksApp2_nEQ=getDistanceEQDist(mapLinksApp2_nEQ,numLinks,mapBridges,0);
            
            % Construct vector of distances, R. Only need one set of
            % params for wEQ and one set of params for nEQ. Open loop or
            % filter does not matter.
            distVec_wEQ=getDistVec(mapLinksApp_wEQ,keyLinks,numBridges,numLinks);
            distVec_nEQ=getDistVec(mapLinksApp_nEQ,keyLinks,numBridges,numLinks);
            
            % Get spectral accelerations
            PGA_wEQ=getPGA(nEQ_wEQ.mag,distVec_wEQ,F,mapLinksApp_wEQ,keyLinks,mapBridges);
            PGA_nEQ=getPGA(nEQ_nEQ.mag,distVec_nEQ,F,mapLinksApp_nEQ,keyLinks,mapBridges);
            
            % Takes the Sa from the attenuation law and gives the corresponding
            % probability of exceedance
            Pex_wEQ=getPex(PGA_wEQ,p,SF,PGAVec);
            Pex_nEQ=getPex(PGA_nEQ,p,SF,PGAVec);

            % Take the probability of exceedance and compute the probability of damage
            % states as a vector
            ProbDS_wEQ=getProbVec(Pex_wEQ);
            ProbDS_nEQ=getProbVec(Pex_nEQ);
            
            % Get the range of the four limit states
            DSlims_wEQ=getDSRange(ProbDS_wEQ);
            DSlims_nEQ=getDSRange(ProbDS_nEQ);
            
            % Assign the parameter based on the probability table (for each
            % ensemble)
            mapLinksApp_wEQ=getMapParams(mapLinksApp_wEQ,numLinks,keyLinks,DSlims_wEQ);
            mapLinksApp2_wEQ=getMapParams(mapLinksApp2_wEQ,numLinks,keyLinks,DSlims_wEQ);
            
            mapLinksApp_nEQ=getMapParams(mapLinksApp_nEQ,numLinks,keyLinks,DSlims_nEQ);
            mapLinksApp2_nEQ=getMapParams(mapLinksApp2_nEQ,numLinks,keyLinks,DSlims_nEQ);
            
        end
        
        % Code to extract qmax values from each link (even those that arent
        % bridges) for each ensemble. 
        
        qmaxArray_wEQ(:,1,i,k+1)=getqmaxMat(mapLinksApp_wEQ,keyLinks);
        qmaxArray2_wEQ(:,1,i,k+1)=getqmaxMat(mapLinksApp2_wEQ,keyLinks);
        
        qmaxArray_nEQ(:,1,i,k+1)=getqmaxMat(mapLinksApp_nEQ,keyLinks);
        qmaxArray2_nEQ(:,1,i,k+1)=getqmaxMat(mapLinksApp2_nEQ,keyLinks);
        
        % DENSITY AND PARAM MODELS.
        
        % Forward progression of densities via CTM
        
        xi_m_wEQ(1:totCells,1,i,k+1)=forwardPredictEns(DT,xi_p_wEQ(:,1,i,k),mapLinksApp_wEQ,...
            mapNodes,totCells,startApp,endApp,1,err_S,err_Q,err_R,Q_S,Q_Q,Q_R,noiseBC);
        
        xi_m2_wEQ(1:totCells,1,i,k+1)=forwardPredictEns(DT,xi_p2_wEQ(:,1,i,k),mapLinksApp2_wEQ,...
            mapNodes,totCells,startApp,endApp,1,err_S,err_Q,err_R,Q_S,Q_Q,Q_R,noiseBC);
        
        xi_m_nEQ(1:totCells,1,i,k+1)=forwardPredictEns(DT,xi_p_nEQ(:,1,i,k),mapLinksApp_nEQ,...
            mapNodes,totCells,startApp,endApp,1,err_S,err_Q,err_R,Q_S,Q_Q,Q_R,noiseBC);
        
        xi_m2_nEQ(1:totCells,1,i,k+1)=forwardPredictEns(DT,xi_p2_nEQ(:,1,i,k),mapLinksApp2_nEQ,...
            mapNodes,totCells,startApp,endApp,1,err_S,err_Q,err_R,Q_S,Q_Q,Q_R,noiseBC);
        
        % Rewrap
        mapLinksApp_wEQ_array{i}=mapLinksApp_wEQ;
        mapLinksApp2_wEQ_array{i}=mapLinksApp2_wEQ;
        mapLinksApp_nEQ_array{i}=mapLinksApp_nEQ;
        mapLinksApp2_nEQ_array{i}=mapLinksApp2_nEQ;
        
    end
    
    % PRIOR MEAN AND ERROR COVARIANCE
    
    % Compute A, and unbiased covariance matrix (Evensen)
    A_wEQ=composeEnsToMat(xi_m_wEQ,totCells,noEn,k+1);
    A_nEQ=composeEnsToMat(xi_m_nEQ,totCells,noEn,k+1);
    
    % Compute ensemble means
    Abar_wEQ=A_wEQ*1/noEn*ones(noEn); 
    Abar_nEQ=A_nEQ*1/noEn*ones(noEn); 
    
    % Put means into desired notation
    xk_m_wEQ(:,1,k+1)=Abar_wEQ(:,1);
    xk_m_nEQ(:,1,k+1)=Abar_nEQ(:,1);
   
    Apr_wEQ=A_wEQ-Abar_wEQ;
    Apr_nEQ=A_nEQ-Abar_nEQ;
    
    Pk_m_wEQ(:,:,k+1)=(1/(noEn-1))*(Apr_wEQ*Apr_wEQ');
    Pk_m_nEQ(:,:,k+1)=(1/(noEn-1))*(Apr_nEQ*Apr_nEQ');
    
    % Open loop
    A2_wEQ=composeEnsToMat(xi_m2_wEQ,totCells,noEn,k+1);
    A2_nEQ=composeEnsToMat(xi_m2_nEQ,totCells,noEn,k+1);
    
    Abar2_wEQ=A2_wEQ*1/noEn*ones(noEn);
    Abar2_nEQ=A2_nEQ*1/noEn*ones(noEn);
    
    xk_m2_wEQ(:,1,k+1)=Abar2_wEQ(:,1);
    xk_m2_nEQ(:,1,k+1)=Abar2_nEQ(:,1);
    
    Apr2_wEQ=A2_wEQ-Abar2_wEQ;
    Apr2_nEQ=A2_nEQ-Abar2_nEQ;
    
    Pk_m2_wEQ(:,:,k+1)=(1/(noEn-1))*(Apr2_wEQ*Apr2_wEQ');
    Pk_m2_nEQ(:,:,k+1)=(1/(noEn-1))*(Apr2_nEQ*Apr2_nEQ');

    % ANALYSIS SCHEME (Evensen)
    
    % Account for damage to the sensor after the EQ occurs
    if k>=keq
        
        % Use the modified measurement matrix
        H=Heq;
        
        % Initialize
        D_wEQ=zeros(length(sensInd),noEn);
        Gam_wEQ=zeros(length(sensInd),noEn);
        
        D_nEQ=zeros(length(sensInd),noEn);
        Gam_nEQ=zeros(length(sensInd),noEn);
             
        for i=1:noEn
            
            meas=y(k+1,sensInd)';
            
            % Get ensemble noise and noisy measurement
            [ensMeas ensMeasNoise]=getEnsMeasNoiseEQ(meas,mapSensors,keySensors,sensInd);
            
            % Add to D and Gam matrices
            D_wEQ(:,i)=ensMeas;
            Gam_wEQ(:,i)=ensMeasNoise;
            D_nEQ(:,i)=ensMeas;
            Gam_nEQ(:,i)=ensMeasNoise;
            
        end

    else
        
        % Initialize
        D_wEQ=zeros(numSensors,noEn);
        Gam_wEQ=zeros(numSensors,noEn);
        
        D_nEQ=zeros(numSensors,noEn);
        Gam_nEQ=zeros(numSensors,noEn);
        
        for i=1:noEn
            
            meas=y(k+1,:)';
            
            % Get ensemble noise and noisy measurement
            [ensMeas ensMeasNoise]=getEnsMeasNoise(meas,mapSensors,keySensors);
            
            % Add to D and Gam matrices
            D_wEQ(:,i)=ensMeas;
            Gam_wEQ(:,i)=ensMeasNoise;
            D_nEQ(:,i)=ensMeas;
            Gam_nEQ(:,i)=ensMeasNoise;
            
        end
        
    end
    
    % Compute Re
    Re_wEQ=(1/(noEn-1))*(Gam_wEQ*Gam_wEQ');
    Re_nEQ=(1/(noEn-1))*(Gam_nEQ*Gam_nEQ');
    
    % Compute D' 
    Dpr_wEQ=D_wEQ-H*A_wEQ;
    Dpr_nEQ=D_nEQ-H*A_nEQ;
    
    % Compute innovation and Kalman gain
    [innov_array_wEQ{k+1},kalGain_array_wEQ{k+1},D_array_wEQ{k+1},Gam_array_wEQ{k+1}]=...
        getUpdateInfo(H,Pk_m_wEQ(:,:,k+1),Re_wEQ,Dpr_wEQ,D_wEQ,Gam_wEQ);
    [innov_array_nEQ{k+1},kalGain_array_nEQ{k+1},D_array_nEQ{k+1},Gam_array_nEQ{k+1}]=...
        getUpdateInfo(H,Pk_m_nEQ(:,:,k+1),Re_nEQ,Dpr_nEQ,D_nEQ,Gam_nEQ);
   
    % Singular value decomposition (SVD)
    
    [U_wEQ,Sig_wEQ,V_wEQ]=svd(H*Apr_wEQ+Gam_wEQ);
    Lam_wEQ=Sig_wEQ*Sig_wEQ';
    
    [U_nEQ,Sig_nEQ,V_nEQ]=svd(H*Apr_nEQ+Gam_nEQ);
    Lam_nEQ=Sig_nEQ*Sig_nEQ';
    
    % Compute update matrix, Aa, from X scheme 
    
    X1_wEQ=Lam_wEQ\U_wEQ';
    X2_wEQ=X1_wEQ*Dpr_wEQ;
    X3_wEQ=U_wEQ*X2_wEQ;
    X4_wEQ=(H*Apr_wEQ)'*X3_wEQ;
    X5_wEQ=eye(noEn)+X4_wEQ;
    
    Aa_wEQ=A_wEQ*X5_wEQ;
    
    X1_nEQ=Lam_nEQ\U_nEQ';
    X2_nEQ=X1_nEQ*Dpr_nEQ;
    X3_nEQ=U_nEQ*X2_nEQ;
    X4_nEQ=(H*Apr_nEQ)'*X3_nEQ;
    X5_nEQ=eye(noEn)+X4_nEQ;
    
    Aa_nEQ=A_nEQ*X5_nEQ;
    
    % Add X5 and Aa to their respective arrays (for smoother)
    X5_array_wEQ{k+1}=X5_wEQ;
    Af_array_wEQ{k+1}=Aa_wEQ;
    
    X5_array_nEQ{k+1}=X5_nEQ;
    Af_array_nEQ{k+1}=Aa_nEQ;
    
    % Do not do analysis scheme!
    Aa2_wEQ=A2_wEQ;
    Aa2_nEQ=A2_nEQ;
    
    % Decompose Aa into ensemble vectors
    
    for i=1:noEn
        
        xi_p_wEQ(:,1,i,k+1)=decompMatToEns(Aa_wEQ,i);
        xi_p2_wEQ(:,1,i,k+1)=decompMatToEns(Aa2_wEQ,i);
        
        xi_p_nEQ(:,1,i,k+1)=decompMatToEns(Aa_nEQ,i);
        xi_p2_nEQ(:,1,i,k+1)=decompMatToEns(Aa2_nEQ,i);
   
    end

    % Compute the unbiased covariance matrix for updated ensembles
    
    Aabar_wEQ=Aa_wEQ*1/noEn*ones(noEn);
    Aabar_nEQ=Aa_nEQ*1/noEn*ones(noEn);
    
    % Put means into desired notation
    xk_p_wEQ(:,1,k+1)=Aabar_wEQ(:,1);
    Aapr_wEQ=Aa_wEQ-Aabar_wEQ;
    Pk_p_wEQ(:,:,k+1)=(1/(noEn-1))*(Aapr_wEQ*Aapr_wEQ');
    
    xk_p_nEQ(:,1,k+1)=Aabar_nEQ(:,1);
    Aapr_nEQ=Aa_nEQ-Aabar_nEQ;
    Pk_p_nEQ(:,:,k+1)=(1/(noEn-1))*(Aapr_nEQ*Aapr_nEQ');
    
    Aabar2_wEQ=Aa2_wEQ*1/noEn*ones(noEn);
    Aabar2_nEQ=Aa2_nEQ*1/noEn*ones(noEn);
    
    % Put means into desired notation
    xk_p2_wEQ(:,1,k+1)=Abar2_wEQ(:,1);
    Aapr2_wEQ=Aa2_wEQ-Aabar2_wEQ;
    Pk_p2_wEQ(:,:,k+1)=(1/(noEn-1))*(Aapr2_wEQ*Aapr2_wEQ');
    
    xk_p2_nEQ(:,1,k+1)=Abar2_nEQ(:,1);
    Aapr2_nEQ=Aa2_nEQ-Aabar2_nEQ;
    Pk_p2_nEQ(:,:,k+1)=(1/(noEn-1))*(Aapr2_nEQ*Aapr2_nEQ');
       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert to matrix form of density means

% WITH EQ
xMatm_wEQ=conv2Mat(xk_m_wEQ); % prior
xMatp_wEQ=conv2Mat(xk_p_wEQ); % posterior
xMatp2_wEQ=conv2Mat(xk_p2_wEQ); % open loop

% WITH NO EQ
xMatm_nEQ=conv2Mat(xk_m_nEQ); % prior
xMatp_nEQ=conv2Mat(xk_p_nEQ); % posterior
xMatp2_nEQ=conv2Mat(xk_p2_nEQ); % open loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc