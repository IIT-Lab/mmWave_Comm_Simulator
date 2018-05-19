%h5gPBCHDMRSIndices PBCH DM-RS resource element indices
% 
%   TS 38.211 section 7.4.3.1.
%  
%   Example:
%   % Generate the 144 resource element indices associated with the PBCH
%   % DM-RS symbol within a single SS/PBCH block.
%
%   indices = h5gPBCHIndices(struct('NCellID',10));
% 
%   See also h5gPBCHDMRS, h5gPBCHIndices, h5gPSSIndices, h5gSSSIndices.

%   Copyright 2018 The MathWorks, Inc.

function indices = h5gPBCHDMRSIndices(gnb)

    vshift = mod(gnb.NCellID,4);
    indices = 1+vshift+[240+(0:4:239), (2*240)+(0:4:47), (2*240)+(192:4:239), (3*240)+(0:4:239)]';
    
end