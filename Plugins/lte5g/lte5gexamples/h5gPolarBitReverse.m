function br = h5gPolarBitReverse(b,n)
%h5gPolarBitReverse Bit-wise reverse input value
%
%   BR = h5gPolarBitReverse(B,N) returns the bit-wise reversed-value of B,
%   each represented over N bits.
%
%   % Example:
%   br = h5gPolarBitReverse(21,7)
%   %  returns br=84
%
%   See also h5gPolarConstruct, h5gPolarEncoder, h5gPolarDecoder.

%   Copyright 2017 The MathWorks, Inc.

%#codegen

% No checking: b, n are scalars, 0<=b<=2^n-1
br = comm.internal.convertBit2Int(flipud( ...
        comm.internal.convertInt2Bit(b,n)),n);

end
