% getUpdateInfo: for a specific time, get the innovation and, Kalman gain,
% D matrix, and Gamma matrix 
%
% INPUTS:
% H: measurement matrix
% Pe: prior error covariance matrix
% Re: matrix of ensemble measurement noise
% Dpr: matrix of ensemble innovation
% D: matrix of ensemble measurement
% Gam: matrix of ensemble noise
%
% NOTE: this inputs follow the Evensen (2003) analysis scheme

function [innov,kalGain,DMat,GamMat]=getUpdateInfo(H,Pe,Re,Dpr,D,Gam)

% Compute Kalman gain
kalGain=Pe*H'/(H*Pe*H'+Re);

innov=Dpr;
DMat=D;
GamMat=Gam;

