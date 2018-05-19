%h5gPSSIndices PSS resource element indices
% 
%   TS 38.211 section 7.4.3.1.
% 
%   Example:
%   % Generate the 127 resource element indices associated with the PSS
%   % within a single SS/PBCH block.
%   
%   indices = h5gPSSIndices();
% 
%   See also h5gPSS, h5gSSSIndices, h5gPBCHIndices.

%   Copyright 2018 The MathWorks, Inc.

function indices = h5gPSSIndices()

    indices = 1+(56:182)'; % 1-based (0th symbol)
    
end
