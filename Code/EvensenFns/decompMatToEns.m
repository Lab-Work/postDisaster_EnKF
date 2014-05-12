% decompMatToEns: takes the matrix holding ensemble members and decomposes
% it to its ensemble vectors
%
% INPUTS
% Aa: matrix holding posterior ensemble members
% i: current ensemble

function ensAVec=decompMatToEns(Aa,i)

ensAVec=Aa(:,i);