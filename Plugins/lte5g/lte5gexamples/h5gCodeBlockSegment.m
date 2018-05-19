%   CBS = h5gCodeBlockSegment(BLK,NBG) splits the input data bit vector BLK
%   into a matrix CBS of code block segments (with filler bits and type-24B
%   CRC appended as appropriate) based on the input base graph number NBG,
%   according to the rules of TS 38.212 5.2.2. Code block segmentation
%   occurs in transport blocks (after CRC appending) for LDPC encoded
%   transport channels (DL-SCH, UL-SCH, PCH).
%
%   The segmentation and padding operation ensures that code blocks
%   entering the LDPC coder are no larger than Kcb (8448 when NBG = 1 or
%   3840 when NBG = 2) in length. If the input block length is > Kcb then
%   the input block is split into a couple of smaller code blocks where
%   each individual block also has a type-24B CRC appended to it. The
%   <NULL> filler bits (represented by -1 at the output) are prepended to
%   the first code block so that all blocks in the set have legal lengths.
%   If the input block length is <= Kcb then no segmentation occurs and no
%   CRC is appended but the single output code block may have <NULL> filler
%   bits prepended. The dimensions of output matrix CBS is K-by-C, where K
%   denotes the length of all code blocks and C the number of code blocks.
%
%   Example:
%   % Code block segmentation occurs if the input length is greater than
%   % the maximum code block size Kcb (8448 when base graph number is 1;
%   %  3840 when base graph number is 2).
%   
%   cbs1 = h5gCodeBlockSegment(ones(4000,1),1); %  No segmentation
%   cbs2 = h5gCodeBlockSegment(ones(4000,1),2); %  With segmentation
%   size(cbs1)
%   size(cbs2)
%
%   % The above example returns:
%   % ans =
%   %         4224           1
%   %
%   %
%   % ans =
%   %
%   %         2080           2
%
%   See also h5gCodeBlockDesegment, h5gLDPCEncode, h5gDLSCH,
%   lteCodeBlockSegment.

%   Copyright 2018 The MathWorks, Inc.

function cbs = h5gCodeBlockSegment(blk,nbg)
    
    % Validate base graph number
    if ~(nbg == 1 || nbg == 2)
        error('lte:error','The input base graph category number (%d) should be either 1 or 2.',nbg);
    end
    
    % Obtain parameters for code block segementation
    chsinfo = getCBSParams(length(blk),nbg);
    
    cbs = zeros(chsinfo.K,chsinfo.C,'int8');
    if chsinfo.C == 1
        blk = int8(blk(:));
        cbs = [blk; -1*ones(chsinfo.F,1,'int8')]; 
    else
        for r = 0:chsinfo.C-1
            cb = blk(r*chsinfo.cbz+1:(r+1)*chsinfo.cbz);
            % Code block CRC encoding
            cb = h5gCRCEncode(cb,'24B');
            % Add filler bits
            cbs(:,r+1) = [cb; -1*ones(chsinfo.F,1,'int8')]; 
        end
    end
    
end

function info = getCBSParams(B,nbg)

    % Get the maximum code block size
    if nbg == 1
        Kcb = 8448;
    else
        Kcb = 3840;
    end
    
    % Get number of code blocks and length of CB-CRC coded block
    if B <= Kcb
        C = 1;
        Bd = B;
    else
        L = 24; % Length of the CRC bits attached to each code block
        C = ceil(B/(Kcb-L));
        Bd = B+C*L;
    end
    
    % Obtain and verify the number of bits per code block (excluding CB-CRC
    % bits)
    cbz = B/C;
    if fix(cbz) ~= (cbz)
        error('lte:error','The length of the input data (%d) is not a multiple integer of the number of code blocks (%d)',B,C);
    end
    
    % Get number of bits in each code block (excluding Filler bits)
    Kd = Bd/C;
    
    % Find the minimum value of Z in all sets of lifting sizes in 38.212
    % Table 5.3.2-1, denoted as Zc, such that Kb*Zc>=Kd
    if nbg == 1
        Kb = 22;
    else
        if B > 640
            Kb = 10;
        elseif B > 560
            Kb = 9;
        elseif B > 192
            Kb = 8;
        else
            Kb = 6;
        end
    end
    Zlist = [2:16 18:2:32 36:4:64 72:8:112 120 128 144 160:16:256 288 320 352 384];
    for Zc = Zlist
        if Kb*Zc >= Kd
            break;
        end
    end
    
    % Get number of bits (including Filler bits) to be encoded by LDPC
    % encoder
    if nbg == 1
        K = 22*Zc;
    else
        K = 10*Zc;
    end
    
    info.C = C;         % Number of code block segements
    info.cbz = cbz;     % Number of bits in each code block (excluding CB-CRC bits and Filler bits)
    info.K = K;         % Number of bits in each code block (including Filler bits)
    info.F = K-Kd;      % Number of fillers bits per code block
    info.ZC = Zc;       % Selected lifting size
    
end