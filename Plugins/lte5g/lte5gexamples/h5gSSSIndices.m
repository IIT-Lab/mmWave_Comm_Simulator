%h5gSSSIndices SSS resource element indices
% 
%   TS 38.211 section 7.4.3.1.
% 
%   Example:
%   % Generate the 127 resource element indices associated with the SSS
%   % within a single SS/PBCH block.
%   
%   indices = h5gSSSIndices();
% 
%   See also h5gSSS, h5gPSSIndices, h5gPBCHIndices.

%   Copyright 2018 The MathWorks, Inc.

function indices = h5gSSSIndices()

    indices = 1+(2*240)+(56:182)'; % 1-based (2th symbol out of 0,1,2,3)
    
end
