% initialCov: gives the initial covariance matrix P0|0 taking into account
% the exact solution, the initial guess, and the noise distribution
%
% INPUTS
% x0acc: vector of true solution initial densities
% x0guess: vector of model solution initial densities
% toler: scale factor
% totCells: integer for number of cells in system

function P0=initialCov(x0acc,x0guess,toler,totCells)

% Initialize
P0=zeros(totCells);

% Populate the covariance matrix 

for i=1:totCells
    sig=(x0acc(i)-x0guess(i))/toler;
    var=sig^2;
    P0(i,i)=var;
end