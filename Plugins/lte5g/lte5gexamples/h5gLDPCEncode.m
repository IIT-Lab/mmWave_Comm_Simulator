%h5gLDPCEncode 5G LDPC encode
%   OUT = h5gLDPCEncode(IN,NBG) returns the result of LDPC encoding in the
%   matrix OUT, based on input data matrix IN and base graph number (1 or
%   2), according to the rules of 38.212 5.3.2. Each column in IN and OUT
%   represents a code block segement before and after LPDC encoding,
%   respectively. The number of columns in IN or OUT equal to the number of
%   code block segments. The encoding procedure includes replacing filler
%   bits (represented by -1) with 0 bits, encoding to generate the full
%   LDPC codeword, replacing the filler bit locations in codeword with -1
%   and systematic bits puncturing.
%
%   Example:
%   % The example below shows a block of 5000 bits is sgemented into two  
%   % code blocks of length 2560, then the LDPC encoder process the two 
%   % code blocks to obtain two encoded code blocks of length 12800.
%
%   blk = ones(5000,1);
%   nbg = 2;
%   cbs = h5gCodeBlockSegment(blk,nbg);
%   codedcbs = h5gLDPCEncode(cbs,nbg);
%   size(cbs)
%   size(codedcbs)
%
%   See also h5gCodeBlockSegment, h5gLDPCDecode, h5gLDPCParityCheckMatrix.

%   Copyright 2017-2018 The MathWorks, Inc.

function out = h5gLDPCEncode(in,nbg)

    % Empty in, empty out
    if isempty(in)
        out = [];
        return;
    end
    
    % Validate base graph number
    if ~(nbg == 1 || nbg == 2)
        error('lte:error','The input base graph category number (%d) should be either 1 or 2.',nbg);
    end
    
    % Obtain input/output size
    K = size(in, 1);
    C = size(in, 2);
    if nbg == 1
        Zc = K/22;
        N = 66*Zc;
    else
        Zc = K/10;
        N = 50*Zc;
    end

    % Validate the input data size
    if fix(Zc) ~= Zc
        if nbg == 1
            nsys = 22;
        else
            nsys = 10;
        end
        error('lte:error','The number of rows in the input data matrix (%d) should be a multiple integer of %d (for base graph %d).', K, nsys, nbg);
    end
    
    % Encode each code block
    out = zeros(N, C, 'int8');
    for r = 1:C
        out(:,r) = singleCBEncode(in(:,r),Zc,nbg);
    end
    
end

% Encode a single code block
function out = singleCBEncode(in,zc,nbg)

    % Initialize persistent parameters
    persistent encoders;
    persistent BGvec;
    persistent Zvec;
    if (isempty(Zvec))
        encoders = cell(2, 51);
        BGvec    = [1 2]; % All possible base graph numbers
        Zvec     = [2:16 18:2:32 36:4:64 72:8:112 120 128 144 160:16:256 288 320 352 384]; % All possible lifting sizes
    end
    
    in = int8(in);
    in = in(:);
    
    % Create encoder on demand if it does not exist
    bgIdx = nbg == BGvec;
    ZIdx = zc == Zvec;
    encoder = encoders{bgIdx, ZIdx};
    if isempty(encoder)
        pcm = h5gLDPCParityCheckMatrix(nbg, zc);
        encoder = comm.LDPCEncoder(sparse(pcm));
        encoders{bgIdx, ZIdx} = encoder;
    end
    
    % Find filler bits (shortening bits) and replace them with 0 bits
    locs = find(in==-1);
    in(locs) = 0;
    
    % LDPC encode
    out = encoder(in);
    
    % Put filler bits back to the decoded bits
    out(locs) = -1;
    
    % Puncture 2*Zc systematic bits
    out(1:2*zc) = [];
    

end

