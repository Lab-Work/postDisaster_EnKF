% inpEQ: Create the earthquake object by specifying the magnitude

function newEQ=inpEQ

prompt={'Magnitude'};

% Input bridges
inp=inputdlg(prompt,'Define EQ',1,{'0'});

newEQ=EQ;

newEQ.mag=str2num(inp{1});