% initialize3DArray: creates 3D arrays
%
% INPUTS
% n: number of cells
% K: number of time steps

function [xk_m xk_p Pk_m Pk_p]=initialize3DArray(n,K)
    
Pk_m=zeros(n,n,K);
xk_m=zeros(n,1,K);
xk_p=zeros(n,1,K);
Pk_p=zeros(n,n,K);