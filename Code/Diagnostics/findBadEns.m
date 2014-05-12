% findBadEns: finds ensembles that are out of range and returns in an
% array that tells the time step, ensemble number, and cell of the bad
% ensemble
%
% INPUTS:
% xi_ens: a 4D array of ensembles
% start: start time step
% fin: end time step
%
% NOTE: The jam density is hardcoded here as 250

function badEnsArray=findBadEns(xi_ens,start,fin)

% Initialize
totCells=size(xi_ens,1);
noEn=size(xi_ens,3);
badEnsArray=[];

% Counter index
count=1;

% Jam density
rhoj=250;

for k=start:fin
    for i=1:noEn
        for c=1:totCells   
            % Look for negative ensembles or ensembles above jam density
            if (xi_ens(c,1,i,k)<0 || xi_ens(c,1,i,k)>rhoj)
                
                badEnsArray(count,1)=k;
                badEnsArray(count,2)=i;
                badEnsArray(count,3)=c;
                
                if xi_ens(c,1,i,k)>rhoj
                    % If above jam density
                    badEnsArray(count,4)=1;   
                else
                    % If negative
                    badEnsArray(count,4)=0;
                end
               
                count=count+1;
                
            end  
        end 
    end
end