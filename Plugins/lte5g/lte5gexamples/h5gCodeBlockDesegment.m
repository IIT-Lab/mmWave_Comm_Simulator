%h5gCodeBlockDesegment Code block desegmentation and CRC decoding
%   [BLK,ERR] = h5gCodeBlockDesegment(...) performs the inverse of the code
%   block segmentation and CRC appending (see <a href="matlab: help('h5gCodeBlockSegment')">h5gCodeBlockSegment</a>). It 
%   concatenates the input code block segments into a single output data
%   block BLK, removing any filler and type-24B CRC bits that may be
%   present in the process. The results of code block CRC decoding (if
%   applicable) are available in vector ERR.
%
%   [BLK,ERR] = h5gCodeBlockDesegment(CBS,NBG,BLKLEN) concatenates the
%   input code blocks contained in matrix CBS into an output vector BLK of
%   length BLKLEN, based on a selected base graph number NBG. BLKLEN is
%   also used to validate the dimensions of the data in CBS and to
%   calculate the amount of filler to be removed. The dimensions of CBS is
%   K-by-C, where K denotes the length of all code blocks and C the number
%   of code blocks. If C > 1 then each code block is assumed to have a
%   type-24B CRC attached. This CRC is decoded and stripped from each code
%   block prior to output concatenation and the CRC error result is placed
%   in the associated element of vector ERR (the length of ERR is equal to
%   the number of code blocks). If C = 1 then no CRC decoding is performed
%   and ERR will be empty. In all cases the number of filler bits removed
%   from the end of a code block is calculated from BLKLEN; if BLKLEN=0
%   then no filler bits are stripped.
%   
%   [BLK,ERR] = h5gCodeBlockDesegment(CBS,NBG) is similar to the above
%   except that no filler bits are stripped from the output. The detection
%   and processing of type-24B CRC is carried out as above.
%   
%   Example:
%   % Code block segmentation occurs if the input length is greater than
%   % 8448 (base graph number is 1). The below example show how 
%   % segmentation and desegmentation happens. The input data of length 
%   % 10000 gets segmented into two code blocks of length 5280, i.e., 
%   % matrix cbs of size 5280-by-2. This undergoes filler bits, CRC 
%   % removal and concatenation resulting in blk of length 10000.
%   
%   nbg = 1;
%   blklen = 10000;
%   cbs = h5gCodeBlockSegment(ones(blklen,1),nbg); 
%   [blk,err] = h5gCodeBlockDesegment(cbs,nbg,blklen);
%   blkSize = size(blk)
%   err
%
%   See also h5gCodeBlockSegment, h5gLDPCDecode, h5gCRCDecode,
%   h5gDLSCHDecode, lteCodeBlockDesegment.

%   Copyright 2018 The MathWorks, Inc.

function [blk, err] = h5gCodeBlockDesegment(cbs,nbg,varargin)

    % Validate the number of inputs
    narginchk(2, 3);
    
    % Validate base graph number
    if ~(nbg == 1 || nbg == 2)
        error('lte:error','The input base graph category number (%d) should be either 1 or 2.',nbg);
    end
    
    % Get number of code blocks
    C = size(cbs, 2);
    
    if nargin == 2 || (nargin == 3 && varargin{1} == 0)
        % no filler bits are stripped from a code block
        chsinfo.F = 0;
    else
        blklen = varargin{1};
        chsinfo = getCBSParams(blklen,nbg);
        
        % Validation dimensions of cbs if there is input for block length
        if C ~= chsinfo.C
            error('lte:error','Number of columns of the input code blocks matrix (%d) conflicts with the one deduced from the input block length (%d).',C,chsinfo.C);
        end
        if size(cbs,1) ~= chsinfo.K
            error('lte:error','Number of rows of the input code blocks matrix (%d) conflicts with the one deduced from the input block length (%d).',size(cbs,1),chsinfo.K);
        end
    end
    
    % Remove filler bits
    cbs(end-chsinfo.F+1:end,:) = [];
    
    if C == 1 % No CB-CRC decoding
        blk = int8(cbs);
        err = [];
        return;
    end

    err = zeros(1,C,'int8');
    blk = zeros(0,1,'int8');
    for r = 1:C
        % CB-CRC decoding
        [cb, err(r)] = h5gCRCDecode(cbs(:,r),'24B');
        % Concatenate decoded code blocks
        blk = [blk; cb]; %#ok<AGROW>
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
    info.K = K;         % Number of bits in each code block (including Filler bits)
    info.F = K-Kd;      % Number of fillers bits per code block
    
end