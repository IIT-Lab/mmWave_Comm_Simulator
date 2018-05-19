%h5gPBCH Physical broadcast channel
%   
%   TS 38.211 Section 7.3.3.
%
%   Example:
%   % Generate the PBCH symbols (QPSK) for the first SS block in a burst
%   (i_SSB=0) from random bits representing encoded BCH bits:
%
%   gnb.NCellID = 17;
%   i_SSB = 0;
%   v = mod(i_SSB,4); % assuming L_max = 4
%   E = 864; % PBCH bit capacity, TS 38.212 Section 7.1.5
%   cw = randi([0 1],E,1);
%
%   sym = h5gPBCH(gnb,v,cw);
%
%   See also h5gPBCHIndices, h5gPBCHDMRS.

%   Copyright 2018 The MathWorks, Inc.

function sym = h5gPBCH(gnb,v,cw)
    
    cinit = gnb.NCellID;
    M_bit = length(cw);
    c = ltePRBS(cinit,[v*M_bit M_bit],'binary');
    scrambled = xor(cw,c);
    
    sym = lteSymbolModulate(scrambled,'QPSK');

end
