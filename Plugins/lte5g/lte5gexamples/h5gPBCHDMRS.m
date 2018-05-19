%h5gPBCHDMRS PBCH demodulation reference signal
% 
%   TS 38.211 section 7.4.1.4.
% 
%   % Example:
%   % Generate the sequence of 144 PBCH DM-RS symbols associated with the
%   % third SSB block in the 2nd half frame of a frame when the SSB
%   % periodicity is 5ms:
%
%   gnb.NCellID = 10;
%   i_SSB = 2; 
%   n_hf = 1;
%   ibar_SSB = (4 * i_SSB) + n_hf;
%
%   dmrs = h5gPBCHDMRS(gnb,ibar_SSB);
% 
%   See also h5gPBCHDMRSIndices, h5gPBCH.

%   Copyright 2018 The MathWorks, Inc.

function dmrs = h5gPBCHDMRS(gnb,ibar_SSB)

    ncellid = gnb.NCellID;
    
    cinit = 2^11*(ibar_SSB+1)*(fix(ncellid/4)+1) + 2^6*(ibar_SSB+1) + mod(ncellid,4);
    prbs = reshape(ltePRBS(cinit,2*144,'signed'),2,144)';
    
    dmrs = (1/sqrt(2))*complex(prbs(:,1),prbs(:,2));
    
end
