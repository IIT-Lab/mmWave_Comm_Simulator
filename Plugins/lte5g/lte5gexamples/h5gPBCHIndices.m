%h5gPBCHIndices PBCH resource element indices
% 
%   TS 38.211 section 7.4.3.1.
% 
%   Example:
%   % Generate the 432 resource element indices associated with the PBCH
%   % symbols within a single SS/PBCH block.
%   
%   indices = h5gPBCHIndices(struct('NCellID',10));
% 
%   See also h5gPBCH, h5gPBCHDMRSIndices, h5gPSSIndices, h5gSSSIndices.

%   Copyright 2018 The MathWorks, Inc.

function [indices,info] = h5gPBCHIndices(gnb)
    
    vshift = mod(gnb.NCellID,4);
    tblock = 1+[1 0 0 0;2 2 1 1;3 3 3 2];
    indices = bsxfun(@plus,[240+(0:4:239), 480+(0:4:47), 480+(192:4:239), 720+(0:4:239)], tblock(:,vshift+1));
    indices = indices(:);
    info = struct('G',numel(indices)*2,'Gd',numel(indices));
    
end
