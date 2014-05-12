% compCap: for a uniform draw from the P(DS) distribution (PMF), return the
% damage state and traffic capacity
%
% INPUTS
% DSlims: DSlims: vector designating ranges of different damage states
% draw: a randomly drawn integer from rand function
% 
% NOTE: This is currently hardcoded for four damage states

function [ratio,DS]=compCap(DSlims,draw)

if DSlims(1)<=draw && draw<DSlims(2)
    ratio=0;
    DS='total';
elseif DSlims(2)<=draw && draw<DSlims(3)
    ratio=.5;
    DS='high';
elseif DSlims(3)<=draw && draw<DSlims(4)
    ratio=.75;
    DS='medium';
elseif DSlims(4)<=draw && draw<DSlims(5)
    ratio=1;
    DS='insignifcant';
else
end