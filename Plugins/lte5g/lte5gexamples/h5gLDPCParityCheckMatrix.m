%h5gLDPCParityCheckMatrix Parity check matrix of 5G LDPC
%   H = h5gLDPCParityCheckMatrix(BGNUM, Z) returns the 5G LDPC parity check
%   matrix H derived from a input base graph number BGNUM and a lifting
%   size Z, according to TS 38.212 5.3.2.
%
%   Example:
%   % Generate the parity check matrix using the selected base graph and
%   % lifting size.
%
%   BGnum = 1;  % Selected base graph 
%   Z = 144;    % Selected lifting size
%   H = h5gLDPCParityCheckMatrix(BGnum, Z);
%   size(H)
%
%   % The above example returns:
%   % 6624        9792
%
%   See also h5gLDPCEncode, h5gLDPCDecode.

%   Copyright 2017-2018 The MathWorks, Inc.

function H = h5gLDPCParityCheckMatrix(BGnum, Z)

    ZSets={[2  4  8  16  32  64 128 256],...
           [3  6 12  24  48  96 192 384],...
           [5 10 20  40  80 160 320],...
           [7 14 28  56 112 224],...
           [9 18 36  72 144 288],...
           [11 22 44  88 176 352],...
           [13 26 52 104 208],...
           [15 30 60 120 240]};

    if ~(BGnum==1 || BGnum==2)
        error(['The input base graph category number (' num2str(BGnum) ') should be either 1 or 2.']);
    end

    % Get lifting set number
    found = false;
    for SetIdx = 1:8 % LDPC lifting size set index
        if any(Z==ZSets{SetIdx})
            found = true;
            break;
        end
    end
    if found == false
        error(['Lifting size set not found for the input lifting size: ' num2str(Z)]);
    end

    % Get the matrix with BG number 'BGnum' and set number 'SetIdx'.
    % The element of matrix V in the following is H_BG(i,j)*V(i,j), where
    % H_BG(i,j) and V(i,j) are defined in TS 38.212 5.3.2; if V(i,j) is not
    % defined in Table 5.3.2-2 or Table 5.3.2-3, the element are -1.
    bgs = load('5GLDPCBaseGraph');
    V = bgs.(['BG' num2str(BGnum) 'S' num2str(SetIdx)]);

    % Get shift values matrix
    % The element of matrix P in the following is P(i,j) in TS 38.212 5.3.2
    % when V(i,j) are defined in Table 5.3.2-2 or Table 5.3.2-3, if not
    % defined, the element are -1.
    P = shiftValCalculation(V, Z);
    
    % Lifting: apply the shift values P for each ZxZ sub-matrix
    % When the element of P is -1, the lifted sub-matrix is zero matrix.
    [numRows, numCols] = size(P);
    H = zeros([numRows numCols]*Z);
    for i = 1:numRows
        for j = 1:numCols
            H(((i-1)*Z+1):i*Z, ((j-1)*Z+1):j*Z) = ...
                BuildSubBlock(Z, P(i,j));
        end
    end
    
end

% Calculate shift values from V for a lifting size Z
function P = shiftValCalculation(V, Z)
    P = zeros(size(V));
    for i = 1:size(V, 1)
        for j = 1:size(V, 2)
            if V(i,j) == -1
                P(i,j) = -1;
            else
                P(i,j) = mod(V(i,j),Z);
            end
        end
    end
end

function M = BuildSubBlock(s, shift)
    if shift == -1
        M = zeros(s);
    else
        M = eye(s);
        if shift > 0
            shift = shift - 1;
            M(:) = M(:,[end-shift:end 1:(end-shift-1)]);
        end
    end
end