%h5gDLSCH Downlink shared channel
%   [CWOUT,CHINFO] = h5gDLSCH(ENB,CHS,OUTLEN,TRBLKIN) applies the complete
%   Downlink Shared Channel (DL-SCH) transport channel coding chain to the
%   input data TRBLKIN and returns the codewords in CWOUT. The processing
%   combines aspects of 3GPP 5G New Radio (NR) working agreements and the
%   existing LTE standard. The encoding process includes type-24A (or
%   type-16) CRC calculation, code block segmentation and type-24B CRC
%   attachment if any, channel encoding (LDPC encoding for 5G, turbo
%   encoding for LTE), rate matching with Redundancy Version (RV) and code
%   block concatenation. The function implements both 5G and LTE DL-SCH
%   encoding processes, the 5G DL-SCH follows TS 38.212 section 5 and the
%   LTE DL-SCH follows TS 36.212 section 5.1. Additional information about
%   the encoding process is returned in the fields of structure CHINFO. The
%   function is capable of processing both a single transport block or
%   pairs of blocks (contained in a cell array) for the case of spatial
%   multiplexing schemes transmitting two codewords. The type of the return
%   variable CWOUT will be the same as input TRBLKIN i.e. if TRBLK is a
%   cell array containing one or two transport blocks then CWOUT will
%   return a cell array of one or two codewords, and if TRBLKIN is a vector
%   of information bits then CWOUT will also return a vector. If encoding a
%   pair of transport blocks then pairs of modulation schemes and RV
%   indicators are required to be defined in the associated parameter
%   fields below. The function does not support code block group
%   (re-)transmissions defined in TS 38.212 5.4.2.1 and TS 38.213 9.1.1.
%   
%   ENB is an input parameter structure containing the fields below, it is
%   required only when CodingType in CHS indicates the LTE DL-SCH encoding
%   chain (CodingType = 'Turbo'). In case of 5G DL-SCH encoding chain
%   (CodingType = 'LDPC'), it is legal for this parameter to be an empty
%   structure.
%   Only required if 'NSoftbits' is provided in CHS:
%      DuplexMode - Optional. Duplex mode ('FDD'(default),'TDD')
%   Only required for 'TDD' duplex mode:
%      TDDConfig  - Optional. Uplink/Downlink Configuration (0...6) 
%                   (default 0)
%   Only required for 'TxDiversity' transmission scheme (as defined in the 
%   parameter below):
%      CellRefP   - Number of cell-specific reference signal antenna ports
%                   (1,2,4).
%
%   CHS is an input parameter structure defining aspects of the Physical
%   Downlink Shared Channel (PDSCH) onto which the codeword(s) will be
%   mapped, and the DL-SCH soft buffer size and redundancy version(s) of
%   the generated codeword(s). The required fields are:
%   CodingType - Optional. Coding type ('LDPC'(default),'Turbo'). It is
%                used to switch between 5G DL-SCH and LTE DL-SCH
%   Modulation - Modulation type character vector or cell array of
%                character vectors (if 2 blocks) associated with each
%                transport block ('QPSK','16QAM','64QAM','256QAM')
%   NLayers    - Total number of transmission layers associated with the
%                transport block(s) (1..8)
%   RV         - Vector of 1 or 2 redundancy version indicators (0,1,2,3)
%
%   When CodingType is 'LDPC', the following additional fields are required:
%   TargetCodeRate - Target coding rate used together with transport block
%                    size to determine the LDPC base graph
%   TBSLBRM        - Optional. Transport block size used to configure the 
%                    soft buffer size for limited buffer rate matching. It 
%                    is assumed no limit is placed on the number of soft
%                    bits, when the field is not provided
%
%   When CodingType is 'Turbo', the following additional fields are
%   optional:
%   TxScheme   - Optional. LTE Transmission schemes, one of:
%                'Port0'       - Single-antenna port, Port 0 (default)
%                'TxDiversity' - Transmit diversity scheme
%                'CDD'         - Large delay CDD scheme
%                'SpatialMux'  - Closed-loop spatial multiplexing scheme
%                'MultiUser'   - Multi-user MIMO scheme
%                'Port5'       - Single-antenna port, Port 5
%                'Port7-8'     - Single-antenna port, port 7 (when 
%                                NLayers=1); Dual layer transmission, port
%                                7 and 8 (when NLayers=2)
%                'Port8'       - Single-antenna port, Port 8
%                'Port7-14'    - Up to 8 layer transmission, ports 7-14
%   NSoftbits  - Optional. Total number of soft buffer bits 
%                (default=0=no buffer limit)
%   
%   OUTLEN is an input vector (one or two elements in length) defining the
%   codeword length(s) that the input transport block(s) should be rate
%   matched to. It represents the PDSCH capacity for the associated
%   codeword and therefore the length(s) of the vector(s) in CWOUT.
%  
%   TRBLKIN is an input parameter containing the transport block(s)
%   information bits to be encoded. It can either be a single vector or a
%   cell array containing one or two vectors. If the latter, then all rate
%   matching calculations assume that the pair will be transmitting on a
%   single PDSCH, distributed across the total number of layers defined in
%   CHS, as per TS 36.211 or TS 38.211. Note that the lowest order
%   information bit of TRBLKIN should be mapped to the most significant bit
%   of the transport block as defined in section 6.1.1 of TS 36.321.
% 
%   CWOUT is the output parameter containing the DL-SCH encoded codewords.
%   It is either a single vector or a cell array containing one or two
%   vectors depending on the type of the input data TRBLKIN.
%
%   The optional output CHINFO is a structure array containing code block
%   segmentation and rate matching related parameters:
%   C       - Total number of code blocks
%   NLayers - Number of layers associated with transport block/codeword
%   Qm      - Bits per symbol variable used in rate matching calculation
%   RV      - RV value associated with one codeword (if RV present 
%             at input)
%
%   When CHS.CodingType is 'LDPC', additional fields are output as follows,
%   the fields also includes the selected base graph used to set up the
%   LDPC codecs.
%   NBG     - Selected base graph number (1 or 2)
%   L       - Length of transport block CRC bits (16 or 24)
%   Kd      - Number of bits in each code block (including code block CRC 
%             bits if any, but excluding Filler bits)
%   K       - Number of bits in each code block (including code block CRC 
%             bits if any, and Filler bits)
%   F       - Number of fillers bits per code block
%   N       - Length of LPDC encoded circular buffer 
%   E       - Rate matching output size for each code blocks
%   Ncb     - Code block soft buffer size
%   ZC      - Selected lifting size
%
%   When CHS.CodingType is 'Turbo', additional fields are output:
%   Km      - Lower code block size (K-)
%   Cm      - Number of code blocks of size Km (C-)
%   Kp      - Upper code block size (K+)
%   Cp      - Number of code blocks of size Kp (C+)
%   F       - Number of filler bits in first block
%   L       - Number of segment CRC bits
%   Bout    - Total number of bits in all segments
%   NL      - Number of layers related variable used in rate matching 
%             calculation
%   NIR     - Number of soft bits associated with transport block
%   
%   Note that if two transport blocks are encoded then CHINFO will be an
%   structure array of two elements, one for each transport block. Note
%   that for LTE, the structure can also be created independently using 
%   the <a href="matlab: help('lteDLSCHInfo')">lteDLSCHInfo</a> function.
%
%   Example:
%   % Generates a 5G DL-SCH codeword with transport block size 28336 and a
%   % taget coding rate of 1/2.
%   
%   enb = struct();
%   chs.Modulation = 'QPSK';
%   chs.NLayers = 1;
%   chs.RV = 0;
%   chs.TargetCodeRate = 0.5;
%   trBlkSize = 28336;
%   codedTrBlkSize = 38880;
% 
%   trBlkData = randi([0,1],trBlkSize,1);
%   codeWord = h5gDLSCH(enb,chs,codedTrBlkSize,trBlkData);
%   codeWord(1:10)
%
%   See also h5gDLSCHInfo, h5gDLSCHDecode, lteDLSCH, lteDLSCHInfo,
%   ltePDSCH.

%   Copyright 2017-2018 The MathWorks, Inc.

function [encbits, chinfo] = h5gDLSCH(enb,dch,outblklen,trblks)
    
    % Validate the number of input arguments
    narginchk(4,4);
    
    % Validate any optional parameters used directly in the MATLAB code
    dch = mwltelibrary('validateLTEParameters',dch,'TxScheme');
    
    % Establish if the input transport data is in a cell array, and if not,
    % place it in a cell array for uniform processing later
    cellout = iscell(trblks);
    if (~cellout)
        trblks = {trblks};
    end
    
    % Validate the input arguments
    outblklen = validateInputs(dch,outblklen,trblks);
    
    % Create the informational output 'chinfo'
    trblklen = cellfun(@length,trblks);
    ldpc = ~isfield(dch,'CodingType') || strcmpi(dch.CodingType,'LDPC'); 
    if ldpc
        chinfo = h5gDLSCHInfo(dch,outblklen,trblklen);
    else
        chinfo = lteDLSCHInfo(enb,dch,trblklen);
    end
    
    % For any empty transport blocks, set the rate matcher output length to
    % zero to ensure that an empty codeword is produced
    outblklen(trblklen==0) = 0;
    
    % Calculate the parameters for encoding each codeword
    chs = createCodewordParameters(dch,chinfo);
    
    % Encode each codeword
    ncw = length(trblks);
    encbits = cell(1,ncw);
    for i = 1:ncw
          encbits{i} = encode(chs(i),outblklen(i),trblks{i});
    end
    
    % If the input transport data was not in a cell array, remove the cell
    % array on the output codeword
    if (~cellout)
        encbits = encbits{1};
    end
    
end

% Encode a single codeword
function cw = encode(chs,outlen,trblkin)
    
    % Optional coding type, default to 5G DL-SCH
    ldpc = strcmpi(chs.CodingType,'LDPC'); 
    if ~ldpc 
        % Transport block CRC attachment
        crced = lteCRCEncode(trblkin,'24A');

        % Code block segmentation and code block CRC attachment
        segmented = lteCodeBlockSegment(crced);

        % Channel coding
        coded = lteTurboEncode(segmented);

        % Rate matching
        cw = lteRateMatchTurbo(coded,outlen,chs.RV,chs);
    else
        % Get transport block CRC type
        A = length(trblkin);
        if A > 3824
            crcPoly = '24A';
        else
            crcPoly = '16';
        end
   
        % Transport block CRC attachment
        crced = h5gCRCEncode(trblkin,crcPoly);
        
        % Code block segmentation and code block CRC attachment
        segmented = h5gCodeBlockSegment(crced,chs.NBG);
        
        % Channel coding
        coded = h5gLDPCEncode(segmented,chs.NBG);

        % Rate matching
        cw = h5gRateMatchLDPC(coded,outlen,chs.RV,chs);
    end
    
end

% Validate the input arguments
function outblklen = validateInputs(dch,outblklen,trblks)

    % Check that relevant parameters are the cell type
    if isstring(dch.Modulation)
        dch.Modulation = arrayfun(@char,dch.Modulation,'UniformOutput',false);
    elseif ~iscell(dch.Modulation)
        dch.Modulation = {dch.Modulation};
    end
    
    % Check that the number of transport blocks presented does not exceed 2
    if length(trblks)>2
        error('lte:error','The number of transport block vectors must be 1 or 2.');
    end
    
    % Check that input cell array (parameter) lengths are the same
    if length(dch.Modulation) < length(trblks)
        error('lte:error','The number of modulation types does not match the number of transport block vectors.');
    end
    if length(outblklen) ~= length(trblks)
        error('lte:error','The number of output codeword lengths (size of OUTLEN) does not match the number of transport block vectors.');
    end
    
    % If the block lengths were parameterized as a cell array then convert to a vector
    if (iscell(outblklen))
        outblklen = cell2mat(outblklen);
    end
    
end

% Create parameters for each codeword, in a structure array 'chs'
function chs = createCodewordParameters(dch,chinfo)
    % The 'modulations' array has modulation strings in the positions of
    % the corresponding number of bits per symbol 'Qm'
    modulations = {'' 'QPSK' '' '16QAM' '' '64QAM' '' '256QAM'};
    modulation = modulations([chinfo.Qm]);
    [chs(1:numel(chinfo)).Modulation] = modulation{:};
    [chs(:).RV] = chinfo.RV;
    [chs(:).NLayers] = chinfo.NLayers;
    if ~isfield(dch,'CodingType')
        dch.CodingType = 'LDPC';
    end
    [chs(:).CodingType] = deal(dch.CodingType);
    if strcmpi(dch.CodingType,'LDPC')
        [chs(:).NBG] = chinfo.NBG;
        if isfield(dch,'TBSLBRM')
            [chs(:).TBSLBRM] = deal(dch.TBSLBRM);
        end
    else
        [chs(:).TxScheme] = deal(dch.TxScheme);
        [chs(:).NIR] = chinfo.NIR;
    end
end