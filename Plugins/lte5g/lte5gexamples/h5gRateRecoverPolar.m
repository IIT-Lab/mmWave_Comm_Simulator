function out = h5gRateRecoverPolar(in,K,N,iBIL)
%h5gRateRecoverPolar Polar rate matching recovery
%
%   OUT = h5gRateRecoverPolar(IN,K,N,IBIL) returns the rate-recovered
%   output, OUT, for an input, IN, of length E. The output, OUT, is of
%   length N and K represents the information block length. IBIL is a 
%   boolean scalar that enables coded bit deinterleaving.
%
%   % Example:
%   N = 2^9;            % Polar encoded block length
%   K = 56;             % Number of information bits
%   E = 864;            % Number of rate matched output bits
%   iBIL = false;       % Deinterleaving of input bits
% 
%   in = randi([0 1],N,1);
%   chIn = h5gRateMatchPolar(in,K,E,iBIL);
%   out = h5gRateRecoverPolar(chIn,K,N,iBIL);
%   isequal(out,in)
%
%   See also h5gRateMatchPolar, h5gPolarDecoder.

%   Copyright 2018 The MathWorks, Inc.

%#codegen

%   Reference:
%   [1] 3GPP TS 38.212, "3rd Generation Partnership Project; Technical 
%   Specification Group Radio Access Network; NR; Multiplexing and channel 
%   coding (Release 15), v15.0.0, 2017-12. Section 5.4.1.

% Channel deinterleaving: Section 5.4.1.3
if iBIL    
    % For uplink only
    inE = iBILDeinterl(in);
else
    % No deinterleaving
    inE = in;
end

% Bit selection: Section 5.4.1.2
E = length(in);
if E >= N 
    % Just the first set output
    outN = inE(1:N);
else
    if K/E <= 7/16
        % puncturing (put at the end)
        outN = zeros(N,1);          % 0s for punctures
        outN(end-E+1:end) = inE;
    else
        % shortening (put at the start)
        outN = 1e20*ones(N,1);      % use a large value for 0s
        outN(1:E) = inE;
    end
end

% Sub-block deinterleaving: Section 5.4.1.1
out = subBlockDeinterleave(outN);

end

%-------------------------------------------------------------------------
function out = subBlockDeinterleave(in)
%
% OUT = subBlockDeinterleave(IN) returns the sub-block deinterleaved input.
% 
%   Reference: TS 38.212, Section 5.4.1.1.

N = length(in);
jn = nr5g.internal.polarSBInterleaveMap(N);
out = zeros(N,1);
out(jn+1) = in;

end

%-------------------------------------------------------------------------
function out = iBILDeinterl(in)
% Triangular deinterleaver
%
%   OUT = iBILDeinterl(IN) performs triangular deinterleaving on the input,
%   IN, and returns the output, OUT.
% 
%   Reference: TS 38.212, Section 5.4.1.3.

% Get T off E
E = length(in);
T = getT(E);

% Create the table with nulls (filled in row-wise)
vTab = zeros(T,T);
k = 0;
for i = 0:T-1
    for j = 0:T-1-i
        if k < E
            vTab(i+1,j+1) = k+1;
        end
        k = k+1;
    end
end

% Write input to buffer column-wise, respecting vTab
v = Inf*ones(T,T);
k = 0;
for j = 0:T-1
    for i = 0:T-1-j
        if k < E && vTab(i+1,j+1) ~= 0
            v(i+1,j+1) = in(k+1);
            k = k+1;
        end
    end
end

% Read output from buffer row-wise
out = zeros(size(in));
k = 0;
for i = 0:T-1
    for j = 0:T-1-i
        if ~isinf(v(i+1,j+1))
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
