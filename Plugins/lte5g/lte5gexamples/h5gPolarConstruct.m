function F = h5gPolarConstruct(K,E,nMax)
%h5gPolarConstruct Polar code construction
%
%   F = h5gPolarConstruct(K,E,NMAX) returns an N-bit vector, F, as the
%   output where K entries in the output would be 0 (information bit
%   locations), and N-K entries in the output would be 1 (frozen bit
%   locations). E is the rate-matched output length and NMAX is the maximum
%   n value (either of 9 or 10). The mother code rate is given by K/N, while
%   the effective code rate after rate-matching is K/E. K, E and NMAX must
%   be all scalars.
%   The function assumes the number of parity check bits, NPC, to be 0.
%
%   % Example: Construct a code with a message length of 48 bits and a 
%   %          rate matched output length of 144.
%   nMax = 9;               % maximum value of n
%   K = 48;                 % message length
%   E = 144;                % code rate
%   F = h5gPolarConstruct(K,E,nMax);
%
%   See also h5gPolarEncoder, h5gRateMatchPolar.

%   Copyright 2017-2018 The MathWorks, Inc.

% References:
%   [1] 3GPP TS 38.212, "3rd Generation Partnership Project; Technical 
%   Specification Group Radio Access Network; NR; Multiplexing and channel 
%   coding (Release 15), v15.0.0, 2017-12. Sections 5.3.1, 5.3.1.2,
%   5.4.1.1.

%#codegen

% K < E
coder.internal.errorIf(K>=E, ...
    'lte5g:h5gPolar:InvalidK_E');

% nMax = [9, 10]
coder.internal.errorIf(~any(nMax==[9 10]), ...
    'lte5g:h5gPolar:InvalidNMax');

% Get n, N, Section 5.3.1
cl2e = ceil(log2(E));
if (E <= (9/8) * 2^(cl2e-1)) && (K/E < 9/16)
    n1 = cl2e-1;
else
    n1 = cl2e;
end

rmin = 1/8;
n2 = ceil(log2(K/rmin));

nMin = 5;
n = max(min([n1 n2 nMax]),nMin);
N = 2^n;

% Get sequence for N, ascending ordered, Section 5.3.1.2
s1024 = nr5g.internal.polarSequence;          % Nmax=1024
idx = (s1024 < N);
qSeq = s1024(idx);                            % 0-based

% Get frozen, information bit indices sets, qF, qI, Section 5.4.1.1
jn = nr5g.internal.polarSBInterleaveMap(N);  % 0-based
qFtmp = [];
if E < N
    if K/E <= 7/16  % puncturing
        for i = 0:(N-E-1)
            qFtmp = [qFtmp; jn(i+1)];
        end
        if E >= 3*N/4
            uLim = ceil(3*N/4-E/2);
            qFtmp = [qFtmp; (0:uLim-1).'];
        else
            uLim = ceil(9*N/16-E/4);
            qFtmp = [qFtmp; (0:uLim-1).'];
        end
        qFtmp = unique(qFtmp);  
    else            % shortening
        for i = E:N-1
            qFtmp = [qFtmp; jn(i+1)];
        end        
    end
end

% Get qI from qFtmp and qSeq
qI = zeros(K,1);            % Assumes npc=0
j = 0;
for i = 1:N
    ind = qSeq(N-i+1);      % flip for most reliable
    if any(ind==qFtmp) 
        continue;
    end
    j = j+1;
    qI(j) = ind;
    if j==K
        break;
    end
end

% Form the frozen bit vector
qF = setdiff(qSeq,qI);    % sorted doesnt matter now
F = zeros(N,1);
F(qF+1) = ones(length(qF),1);

end