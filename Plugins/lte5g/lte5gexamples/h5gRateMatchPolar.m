function out = h5gRateMatchPolar(in,K,E,iBIL)
%h5gRateMatchPolar Polar rate matching
%
%   OUT = h5gRateMatchPolar(IN,K,E,IBIL) returns the rate-matched output,
%   OUT, for a polar encoded input, IN, for an information block length K.
%   The output is of length E. IBIL is a boolean scalar that enables coded
%   bit interleaving. The input IN must be a column vector of N elements,
%   where N is a power of two and greater than K.
%
%   % Example:
%   N = 2^9;            % Polar encoded block length
%   K = 56;             % Number of information bits
%   E = 864;            % Number of rate-matched output bits
%   iBIL = false;       % Interleaving of rate-matched coded bits
% 
%   in = randi([0 1],N,1); 
%   out = h5gRateMatchPolar(in,K,E,iBIL);
%
%   See also h5gRateRecoverPolar, h5gPolarEncoder.

%   Copyright 2018 The MathWorks, Inc.

%#codegen

%   Reference:
%   [1] 3GPP TS 38.212, "3rd Generation Partnership Project; Technical 
%   Specification Group Radio Access Network; NR; Multiplexing and channel 
%   coding (Release 15), v15.0.0, 2017-12. Section 5.4.1.

% Sub-block interleaving: Section 5.4.1.1
y = subBlockInterleave(in);

% Bit selection: Section 5.4.1.2
N = length(in);
outE = zeros(E,1);
if E >= N 
    % Bit repetition
    for k = 0:E-1
        outE(k+1) = y(mod(k,N)+1);
    end
else
    if K/E <= 7/16
        % puncturing (take from the end)
        outE = y(end-E+1:end);
    else
        % shortening (take from the start)
        outE = y(1:E);
    end
end

% Interleaving: Section 5.4.1.3
if iBIL
    % For uplink only
    out = iBILInterl(outE);
else
    % No interleaving
    out = outE;
end

end

%-------------------------------------------------------------------------
function out = subBlockInterleave(in)
%
% OUT = subBlockInterleave(IN) returns the sub-block interleaved output.
% 
%   Reference: TS 38.212, Section 5.4.1.1.
   
N = length(in);
jn = nr5g.internal.polarSBInterleaveMap(N);
out = in(jn+1);

end

%-------------------------------------------------------------------------
function out = iBILInterl(in)
% Triangular interleaver
%
%   OUT = iBILInterl(IN) performs triangular interleaving on the input, IN,
%   writing in the input E elements row-wise and returns the output, OUT,
%   by reading them out column-wise.
% 
%   Reference: TS 38.212, Section 5.4.1.3.

% Get T off E
E = length(in);
T = getT(E);

% Write input to buffer row-wise
v = -1*ones(T,T);   % filler bits
k = 0;
for i = 0:T-1
    for j = 0:T-1-i
        if k < E
            v(i+1,j+1) = in(k+1);
        end
        k = k+1;
    end
end

% Read output from buffer column-wise
out = zeros(size(in));
k = 0;
for j = 0:T-1
    for i = 0:T-1-j
        if v(i+1,j+1) ~= -1
            out(k+1) = v(i+1,j+1);
            k = k+1;
        end
    end
end

end

%--------------------------------------------------------------------
function t = getT(E)
% Use quadratic solution with ceil for >= in expression.

t = ceil((-1+sqrt(1+8*E))/2);

end
