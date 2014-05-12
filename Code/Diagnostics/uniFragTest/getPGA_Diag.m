% getPGA_Diag: takes in magnitude, distance, and site parameters and gives a
% vector of PGA values as dervied using the Campbell 1997 attenuation law
% (DIAGNOSTICS).
%
% INPUTS
% M: an integer for the magnitude
% R: a vector of distances
% F: an integer of the fault type
% S_SR: a vector of soil parameters
% S_HR: a vector of soil parameters

function PGA=getPGA_Diag(M,R,F,S_SR,S_HR)

% Initialize
PGA=zeros(length(R),1);

% Attenuation Relationship by Campbell (1997)
for i=1:length(R)
    
lnPGA(i)=-3.512+0.904*M-1.328*log(sqrt(R(i)^2+(0.149*exp(0.647*M))^2))+...
    (1.125-0.112*log(R(i))-0.0957*M)*F+(0.440-0.171*log(R(i)))*S_SR(i)+...
    (0.405-0.222*log(R(i)))*S_HR(i);

% Convert lnPGA to PGA (units of g)
PGA(i)=exp(lnPGA(i));

end
