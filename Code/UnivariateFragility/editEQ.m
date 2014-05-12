% editEQ: creates a new EQ object by overwriting the old one.
%
% INPUTS
% m: magnitude

function eqObj=editEQ(m)

% New object
eqObj=EQ;

% Assign characteristics
eqObj.mag=m;
