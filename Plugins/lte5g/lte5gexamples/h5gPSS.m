%h5gPSS Primary synchronization signal
% 
%   TS 38.211 section 7.4.2.2.
% 
%   Example:
%   % Generate the 127 primary synchronization signal symbols (BPSK) for a 
%   % given cell ID. The PSS is transmitted in symbol #0 of a SS/PBCH 
%   % block.
%   
%   pss = h5gPSS(struct('NCellID',1));
% 
%   See also h5gPSSIndices, h5gSSS, h5gPBCH.

%   Copyright 2018 The MathWorks, Inc.

function pss = h5gPSS(gnb) 

    dpss = [...
     1  -1  -1   1  -1  -1  -1  -1   1   1  -1  -1  -1   1   1  -1   1  -1   1 ...
    -1  -1   1   1  -1  -1   1   1   1   1   1  -1  -1   1  -1  -1   1  -1   1 ...
    -1  -1  -1   1  -1   1   1   1  -1  -1   1   1  -1   1   1   1  -1   1   1 ...
     1   1   1   1  -1   1   1  -1   1   1  -1  -1   1  -1   1   1  -1  -1  -1 ...
    -1   1  -1  -1  -1   1   1   1   1  -1  -1  -1  -1  -1  -1  -1   1   1   1 ...
    -1  -1  -1   1  -1  -1   1   1   1  -1   1  -1   1   1  -1   1  -1  -1  -1 ...
    -1  -1   1  -1   1  -1   1  -1   1   1   1   1  -1].';

    % Primary synchronization signal
    ncellid = gnb.NCellID;   % 1008 unique physical-layer cell identities 
    n2off = 43*mod(ncellid,3);
    pss = dpss(1+mod(n2off:n2off+126,127));

end
