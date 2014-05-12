% multiRun_Main: main file for executing multiple simulation runs with the
% same initial conditions without parallel computing.

close all; clear all; clc

% Number of simulations
numSims=3;

% Number of estimators to compute BEEQ of
numQuants=3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for simCount=1:numSims
    
    tic
    
    % Run the script
    EQTrafficModel_Main_EQandNoEQ_multiRun
    
    % Compute ratios based on BEEQ formuation
    
    ratioVec=zeros(1,numQuants);
    
    % No EQ, Filter
    ratioVec(1)=computeBEEQ(xMatNoise,xMatp2_nEQ,xMatp_nEQ);
    
    % EQ, open loop
    ratioVec(2)=computeBEEQ(xMatNoise,xMatp2_nEQ,xMatp2_wEQ);
    
    % EQ, filter
    ratioVec(3)=computeBEEQ(xMatNoise,xMatp2_nEQ,xMatp_wEQ);
    
    ratioMat(simCount,:)=ratioVec;
    
    % Store a the time per simulation in a vector
    simTime=toc;
    simtVec(simCount)=simTime;
    
    disp(['Completed sim ' num2str(simCount) ' out of ' num2str(numSims) ', time of sim: ' num2str(simTime)...
        ', total time: ' num2str(sum(simtVec))]);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To compute true BEEQ of MC sims for the estimators

trueBEEQ=zeros(1,numQuants);

for i=1:numQuants
    
    for j=1:numSims
        
        trueBEEQ(i)=trueBEEQ(i)+log10(ratioMat(j,i));
        
    end
    
end

ratioMat
trueBEEQ=(1/numSims)*trueBEEQ;
trueBEEQ=10.^trueBEEQ