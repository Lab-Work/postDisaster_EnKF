% getFragCurves_Diag: takes in the median PGA and dispersion parameters for
% different limit states and creates fragility curves (DIAGNOSTICS)
%
% INPUTS
%
% PGAVec: vector of PGA values for fragility model (defines range)
% med: vector of median PGA values for different limit states
% disp: vector of dispersion values for different limit states

function p=getFragCurves_Diag(PGAVec,med,dsp)

p=zeros(length(med),length(PGAVec));

% Fragility curves 
for i=1:length(med)
        
        p(i,:)=normcdf(log(PGAVec),log(med(i)),dsp(i));
        
end