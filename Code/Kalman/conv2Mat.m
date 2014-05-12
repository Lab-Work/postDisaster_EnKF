% conv2Mat: convert the 3D array into a matrix
%
% INPUTS:
% xk: vector of densities

function xMat=conv2Mat(xk)

numRows=size(xk,3);
numCols=size(xk,1);

% Initialize
xMat=zeros(numRows,numCols);

for i=1:numRows
    
    currRow=xk(:,1,i);
    
    xMat(i,:)=currRow;
    
end
