%h5gLDPCDecode 5G LDPC decode
%   OUT = h5gLDPCDecode(IN,NBG,NLDPCDECITS) performs the inverse of the code
%   LPDC encoding (see <a
%   href="matlab:help('h5gLDPCEncode')">h5gLDPCEncode</a>), it returns the result of LDPC decoding
%   in a matrix OUT, based on input data matrix IN and base graph number (1
%   or 2). Each column in IN and OUT represents a code block segement
%   before and after LPDC decoder, respectively. The number of columns in
%   IN or OUT equal to the number of code block segments. The decoding
%   iterations is terminated if the parity check is satisfied, NLDPCDECITS
%   defines the maximum number of iterations in the decoder, its default
%   value is 25. The decoding procedure includes puncturing bits padding
%   and decoding to get the complete LDPC systematic bits (including the
%   puncturing bits).
%
%   Example:
%   % The example below shows a block of 5000 bits is sgemented into two  
%   % code blocks of length 2560, then the LDPC encoder process the two 
%   % code blocks to obtain two encoded code blocks of length 12800, the
%   % encoded bits are then converted to soft bits to be decoded by the 
%   % LDPC decoder. Finally, the code blocks before encoding is compared 
%   % with decoded code blocks.
%
%   blklen = 5000;
%   nbg = 2;
%   blk = ones(blklen,1);
%   txcbs = h5gCodeBlockSegment(blk,nbg);
%   txcodedcbs = h5gLDPCEncode(txcbs,nbg);
% 
%   % Convert transmitted coded code blocks into soft bits
%   FillerIndices = find(txcodedcbs(:,1) == -1);
%   txcodedcbs(FillerIndices,:) = 0; % Filler bits are treated as 0s when perform encoding
%   rxcodedcbs = double(1-2*txcodedcbs);
% 
%   % Perform decoding with 25 iterations
%   rxcbs = h5gLDPCDecode(rxcodedcbs,nbg,25);
% 
%   F = length(FillerIndices); % Number of filler bits
%   txcbs(end-F+1:end,:) = 0;  % Replace filler bits with 0
% 
%   isequal(rxcbs,txcbs)
%
%   See also h5gCodeBlockSegment, h5gLDPCEncode, h5gLDPCParityCheckMatrix.

%   Copyright 2017-2018 The MathWorks, Inc.

function out = h5gLDPCDecode(in,nbg,nldpcdecits)

    % Validate the input base graph number
    if ~(nbg == 1 || nbg == 2)
        error('lte:error','The input base graph category number (%d) should be either 1 or 2.',nbg);
    end
    
    % Ouput empty is the input data is empty 
    if isempty(in)
        out = [];
        return;
    end
    
    % Obtain input/output size
    N = size(in, 1);
    C = size(in, 2);
    if nbg == 1
        Zc = N/66;
        K = 22*Zc;
    else
        Zc = N/50;
        K = 10*Zc;
    end
    
    % Validate the input data size
    if fix(Zc) ~= Zc
        if nbg == 1
            ncwnodes = 66;
        else
            ncwnodes = 50;
        end
        error('lte:error','The number of rows in the input data matrix (%d) should be a multiple integer of %d (for base graph %d).', N, ncwnodes, nbg);
    end
    
    % Decode each code block
    out = zeros(K, C, 'int8');
    for r = 1:C
        out(:,r) = singleCBDecode(in(:,r),Zc,nbg,nldpcdecits);
    end

end

% Decode a single code block
function out = singleCBDecode(in,zc,nbg,nldpcdecits)

    % Initialize persistent parameters
    persistent decoders;
    persistent BGvec;
    persistent Zvec;
    persistent NDecIts;
    
    if (isempty(Zvec))
        decoders = cell(2, 51, 20);
        BGvec    = [1 2]; % All possible base graph categories
        Zvec     = [2:16 18:2:32 36:4:64 72:8:112 120 128 144 160:16:256 288 320 352 384]; % All possible lifting sizes
        NDecIts  = 5:5:100; % Number of decoding iterations
    end
    
    in = in(:);
    
    % Add punctured 2*Zc systematic bits
    in = [zeros(2*zc,1); in];
    
    % Create decoder if it does not exist
    bgIdx = nbg == BGvec; 
    ZIdx = zc == Zvec; 
    NDecItsIdx = nldpcdecits == NDecIts;
    if ~any(NDecItsIdx)
        error('lte:error','The number of LDPC decoding iterations should be one of 5, 10, 15, 20,..., 100.');
    end
    decoder = decoders{bgIdx, ZIdx, NDecItsIdx};
    if isempty(decoder)
        pcm = h5gLDPCParityCheckMatrix(nbg, zc);
        decoder = comm.LDPCDecoder(sparse(pcm));
        decoder.MaximumIterationCount = NDecIts(NDecItsIdx);
        decoder.IterationTerminationCondition = 'Parity check satisfied';
        decoders{bgIdx, ZIdx, NDecItsIdx} = decoder;
    end
    
    % Decode the soft-bits
    out = int8(decoder(in));
    
end