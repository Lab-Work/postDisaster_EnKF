% initialize4DArray: creates 4D arrays
%
% INPUTS
% n: number of cells
% noEn: number of ensembles
% K: number of time steps

function [xi_m xi_p]=initialize4DArray(n,noEn,K)
    
xi_m=zeros(n,1,noEn,K);
xi_p=zeros(n,1,noEn,K);
