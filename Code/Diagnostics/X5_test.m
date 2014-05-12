% X5_test: Checks each column of an X5 matrix (defined by Evensen) to see
% if it passes the following tests:
% 
% 1. To be unbiased, the sum of each column of X5 should be equal to 1
% 2. X5 should be diagonal dominant in most cases (plot)
% 3. There will be negative off-diagonal terms
%
% INPUTS
% X5: this is X5 matrix, not the array X5_array

function X5_test(X5)

numCols=size(X5,2);

% Compute sum of columns

for i=1:numCols
    
    s=sum(X5(:,i));
    m=max(X5(:,i));
    
    disp(['Column: ' num2str(i)  ', Sum: ' num2str(s) ', Loc of Max: '...
        num2str(find(X5(:,i)==m))]);
    
end

% Plot

figure

imagesc(X5);
colorbar