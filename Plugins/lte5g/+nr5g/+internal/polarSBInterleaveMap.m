function jn = polarSBInterleaveMap(N)
%Subblock interleaving pattern for Polar rate-matching
%
%   out = nr5g.internal.polarSBInterleaveMap(N) returns the sub-block
%   interleaving pattern for length N.
% 
%   See also h5gRateMatchPolar, h5gRateRecoverPolar.

%   Copyright 2018 The MathWorks, Inc.
   
%#codegen

%   Reference:
%   [1] 3GPP TS 38.212, "3rd Generation Partnership Project; Technical 
%   Specification Group Radio Access Network; NR; Multiplexing and channel 
%   coding (Release 15), v15.0.0, 2017-12. Section 5.4.1.1.

% Table 5.4.1.1-1: Sub-block interleaver pattern
pi = [0;1;2;4; 3;5;6;7; 8;16;9;17; 10;18;11;19;
      12;20;13;21; 14;22;15;23; 24;25;26;28; 27;29;30;31];
  
jn = zeros(N,1);
for n = 0:N-1
    i = floor(32*n/N);
    jn(n+1) = pi(i+1)*(N/32)+ mod(n,N/32);
end

end
