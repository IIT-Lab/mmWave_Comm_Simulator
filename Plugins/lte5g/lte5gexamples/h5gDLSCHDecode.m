%h5gDLSCHDecode Downlink shared channel decoding
%   [TRBLKOUT,BLKCRC,STATEOUT] =
%   h5gDLSCHDecode(ENB,CHS,TRBLKLEN,CWIN,STATEIN) returns the information
%   bits TRBLKOUT decoded from the input soft LLR codeword data CWIN. The
%   Downlink Shared Channel (DL-SCH) decoder includes rate recovery,
%   channel decoding (LDPC decoding for 5G, turbo decoding for LTE), code
%   block concatenation and CRC calculations. The function allows 5G and
%   LTE DL-SCH decoding processes, which use the inverse of the operations 
%   described in <a href="matlab:help('h5gDLSCH')">h5gDLSCH</a>. The function also returns the type-24A (or type-16) 
%   transport block CRC decoding result in BLKCRC and the HARQ process
%   decoding state in STATEOUT. The initial HARQ process state can be input
%   via the optional STATEIN parameter. The function is capable of
%   processing both a single codeword or pairs of codewords (contained in a
%   cell array) for the case of spatial multiplexing schemes transmitting
%   two codewords. The type of the return variable TRBLKOUT will be the
%   same as input CWIN, i.e. if CWIN is a cell array containing one or two
%   codewords then TRBLKOUT will return a cell array of one or two
%   transport blocks, and if CWIN is a vector of soft data then TRBLKOUT
%   will also return a vector. If decoding a pair of codewords then pairs
%   of modulation schemes and Redundancy Version (RV) indicators are
%   required to be defined in the associated parameter fields below. The
%   function does not support code block group (re-)transmissions defined
%   in TS 38.212 5.4.2.1 and TS 38.213 9.1.1.
%   
%   ENB is an input parameter structure containing the fields below, it is
%   required only when CodingType in CHS indicates the LTE DL-SCH decoding
%   chain (CodingType = 'Turbo'). In case of 5G DL-SCH decoding chain
%   (CodingType = 'LDPC'), it is legal for this parameter to be an empty
%   structure.
%   Only required if 'NSoftbits' is provided in CHS:
%      DuplexMode - Optional. Duplex mode ('FDD'(default),'TDD')
%   Only required for 'TDD' duplex mode:
%      TDDConfig  - Optional. Uplink/Downlink Configuration (0...6) 
%                   (default 0)
%
%   As the duplex mode will default to FDD if the 'DuplexMode' field is
%   absent, it is legal for this parameter to be an empty structure.
%
%   CHS is an input parameter structure defining aspects of the Physical
%   Downlink Shared Channel (PDSCH) onto which the codeword(s) were be
%   mapped, and the DL-SCH soft buffer size and redundancy version(s) of
%   the received codeword(s). The required fields are:
%   CodingType - Optional. Coding type ('LDPC'(default),'Turbo'). It is
%                used to switch between 5G DL-SCH and LTE DL-SCH
%   Modulation - Modulation type character vector or cell array of
%                character vectors (if 2 transport blocks) associated with 
%                each transport block ('QPSK','16QAM','64QAM','256QAM')
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
%   NLDPCDecIts    - Number of LDPC decoder iteration cycles
%                    (5, 10, 15, 20,...,100) (default 25)
%
%   When CodingType is 'Turbo', the following additional fields are
%   optional:
%   TxScheme   -  Optional. Transmission scheme, one of:
%                 'Port0'       - Single-antenna port, Port 0 (default)
%                 'TxDiversity' - Transmit diversity scheme
%                 'CDD'         - Large delay CDD scheme
%                 'SpatialMux'  - Closed-loop spatial multiplexing scheme
%                 'MultiUser'   - Multi-user MIMO scheme
%                 'Port5'       - Single-antenna port, Port 5
%                 'Port7-8'     - Single-antenna port, port 7 (when 
%                                 NLayers=1); Dual layer transmission, port
%                                 7 and 8 (when NLayers=2)
%                 'Port8'       - Single-antenna port, Port 8
%                 'Port7-14'    - Up to 8 layer transmission, ports 7-14
%   NSoftbits    - Optional. Total number of soft buffer bits 
%                  (default=0=no buffer limit)
%   NTurboDecIts - Optional. Number of turbo decoder iteration cycles 
%                  (1...30) (default 5)
%   
%   TRBLKLEN is an input vector (one or two elements in length) defining
%   the transport block length(s) that the input code block(s) should be
%   rate recovered and decoded to.
%  
%   CWIN is an input parameter containing the floating point soft LLR data
%   of the codeword(s) to be decoded. It can either be a single vector or a
%   cell array containing one or two vectors. If the latter, then all rate
%   matching calculations assume that the pair were transmitting on a
%   single PDSCH, distributed across the total number of layers defined in
%   CHS, as per TS 36.211 or TS 38.211.
% 
%   STATEIN is an optional input structure array (empty or one or two
%   elements) which can input the current decoder buffer state for each
%   transport block in an active HARQ process. If STATEIN is not an empty
%   array and it contains a non-empty field CBSBuffers then this field
%   should contain a cell array of vectors representing the LLR soft buffer
%   states for the set of code blocks at the input to the LDPC or turbo
%   decoder i.e. after explicit rate recovery. The updated buffer states
%   after decoding are returned in the CBSBuffers field in the output
%   parameter STATEOUT. The STATEIN array would normally be generated and
%   recycled from the STATEOUT of previous calls to h5gDLSCHDecode as part
%   of a sequence of HARQ transmissions.
%   
%   TRBLKOUT is the output parameter containing the decoded information
%   bits. It is either a single vector or a cell array containing one
%   or two vectors, depending on the class and dimensionality of CWIN.
%
%   BLKCRC is an output array (one or two elements) containing the result
%   of the type-24A (or type-16) transport block CRC decoding for the
%   transport block(s).
%   
%   STATEOUT, the final output parameter, is a one or two element structure
%   array containing the internal state of each transport block decoder in
%   the following fields:
%   CBSBuffers - Cell array of vectors representing the LLR soft buffer 
%                states for the set of code blocks associated with a single
%                transport block; the buffers are positioned at the input 
%                to the LDPC or turbo decoder i.e. after explicit rate 
%                recovery
%   CBSCRC     - Array of type-24B code block set CRC decoding results
%   BLKCRC     - Type-24A (or type-16) transport block CRC decoding error
%
%   The STATEOUT array would normally be reapplied via the STATEIN variable
%   of subsequent h5gDLSCHDecode function calls as part of a sequence of
%   HARQ retransmissions.
%
%   Example:
%   % Generate and decode 2 transmissions (RV=0 then RV=1) as part of a 
%   % single 5G DL-SCH codeword HARQ process.
%
%   enb = struct();
%   chs.Modulation = 'QPSK';
%   chs.NLayers = 1;
%   chs.TargetCodeRate = 0.5;
%   trBlkSize = 28336;
%   codedTrBlkSize = 38880;
% 
%   trBlkData = randi([0,1],trBlkSize,1);
%   chs.RV = 0; % Create a codeword with RV = 0
%   cw = h5gDLSCH(enb,chs,codedTrBlkSize,trBlkData);
%   cw(cw == 0) = -1;        % Turn logical bits into 'LLR' data
%   % Initialize the decoder states for the first HARQ transmission 
%   decState = [];
%   [rxTrBlk,~,decState] = h5gDLSCHDecode(enb,chs,trBlkSize,cw,decState);
%   % Returned decState contains the decoder buffer state for each 
%   % transport block for an active HARQ process
%   % Create a second retransmitted codeword with RV = 1
%   chs.RV = 1;        
%   cw = h5gDLSCH(enb,chs,codedTrBlkSize,trBlkData);
%   cw(cw == 0) = -1;  % Turn logical bits into 'LLR' data
%   % Previous transmission decoder buffer state, decState, is used
%   % as part of the sequence of active HARQ transmissions
%   rxTrBlk = h5gDLSCHDecode(enb,chs,trBlkSize,cw,decState);
%
%   See also h5gDLSCHInfo, h5gPDSCHDecode, h5gDLSCH, lteDLSCHInfo.

%   Copyright 2017-2018 The MathWorks, Inc.

function [out, err, dstate] = h5gDLSCHDecode(enb,dch,trblklen,sbits,dstate)

    % Validate the number of input arguments
    narginchk(4, 5);
    
    % Create an empty decoder state if it was not an input
    if nargin < 5
        dstate = [];
    end
    
    % Validate any optional parameters used directly in the MATLAB code
    dch = mwltelibrary('validateLTEParameters',dch,'TxScheme','NTurboDecIts');
    
    % Establish if the input codeword is in a cell array (for the case of a
    % single codeword), and if not, place it in a cell array for uniform
    % processing later
    cellout = iscell(sbits);
    if (~cellout)
        sbits = {sbits};
    end
    
    % Validate the input arguments
    trblklen = validateInputs(dch,trblklen,sbits);
    
    % Create the informational output 'chinfo', note that CellRefP only
    % affects chinfo.NLayers for the TxDiversity scheme, and in turn
    % lteRateRecoverTurbo does not use NLayers for the TxDiversity scheme,
    % so we can provide a value for CellRefP here to avoid lteDLSCHDecode
    % requiring it
    if (~isfield(enb,'CellRefP'))
        enb.CellRefP = 1;
    end
    ldpc = ~isfield(dch,'CodingType') || strcmpi(dch.CodingType,'LDPC'); 
    if ldpc
        outblklen = cellfun(@length,sbits);
        chinfo = h5gDLSCHInfo(dch,outblklen,trblklen);
    else
        chinfo = lteDLSCHInfo(enb,dch,trblklen(1:numel(sbits)));
    end
    
    % Calculate the parameters for decoding each codeword
    chs = createCodewordParameters(dch,chinfo);
    
    % Prepare the dstate structure vector
    ncw = length(sbits);
    if isempty(dstate)
        [dstate(1:ncw).CBSBuffers] = deal({});
    elseif ~isstruct(dstate)
         error('lte:error','The decoder state must be a structure array.');
    end
    if length(dstate) < ncw
        [dstate(end+1:ncw).CBSBuffers] = deal({});
    elseif ~isfield(dstate,'CBSBuffers')
        [dstate.CBSBuffers] = deal({});
    end
    
    % Decode each codeword
    out = cell(1,ncw);
    err = ones(1,ncw,'logical');
    for i = 1:ncw
        [out{i},err(i),dstate(i).CBSBuffers,dstate(i).CBSCRC] = decode(chs(i),trblklen(i),sbits{i},dstate(i).CBSBuffers);
        dstate(i).BLKCRC = err(i);
    end
    
    % If the input codeword was not in a cell array (for the case of a
    % single codeword), remove the cell array on the output transport block
    if (~cellout)
        out = out{1};
    end
    
end

% Decode a single codeword
function [out,err,cbsbuffers,segerr] = decode(chs,trblklen,sbits,cbsbuffers)

    % Optional coding type, default to 5G DL-SCH
    ldpc = strcmpi(chs.CodingType,'LDPC'); 
    if ~ldpc
        
        % Rate recovery
        raterecovered = lteRateRecoverTurbo(sbits,trblklen,chs.RV,chs,cbsbuffers);
    
        % Channel decoding
        decoded = lteTurboDecode(raterecovered,chs.NTurboDecIts);
        
        % Code block desegmentation and code block CRC decoding
        [desegmented,segerr] = lteCodeBlockDesegment(decoded,trblklen+24);
        segerr = (segerr~=0);

        % Transport block CRC decoding
        [out,err] = lteCRCDecode(desegmented,'24A');
        
    else
        % To maintain the same as that in the turbo decoding process
        sbits = -1*sbits;
        
        % Base graph selection
        if trblklen > 3824
            crcPoly = '24A';
        else
            crcPoly = '16';
        end
        
        % Rate recovery
        raterecovered = h5gRateRecoverLDPC(sbits,trblklen,chs.RV,chs,cbsbuffers);

        % Channel decoding
        decoded = h5gLDPCDecode(raterecovered,chs.NBG,chs.NLDPCDecIts);
        
        % Code block desegmentation and code block CRC decoding
        [desegmented,segerr] = h5gCodeBlockDesegment(decoded,chs.NBG,trblklen+str2double(crcPoly(1:2)));
        segerr = (segerr~=0);
        
        % Transport block CRC decoding
        [out,err] = h5gCRCDecode(desegmented,crcPoly);
    end
    
    if (isempty(out))
        % treat decoding to a TBS of zero as a non-existent rather than
        % empty transport block, so set the CRC to a pass
        err = false;
    else
        err = (err~=0);
    end
    
    % For an empty data input, the code block soft buffer state should be
    % empty, otherwise it is the output of the rate recovery
    % rate recovery
    if (isempty(sbits))
        cbsbuffers = cell(1,0);
    else
        cbsbuffers = raterecovered;
    end
    
end

% Validate the input arguments
function trblklen = validateInputs(dch,trblklen,sbits)

    % Check that relevant parameters are the cell type
    if ~iscell(dch.Modulation)
        dch.Modulation = {dch.Modulation};
    end
    
    % Check that the number of codewords presented does not exceed 2
    if length(sbits)>2
        error('lte:error','The number of codeword vectors must be 1 or 2.');
    end
    
    % Check that input cell array (parameter) lengths are the same
    if length(dch.Modulation) < length(sbits)
        error('lte:error','The number of modulation types does not match the number of codeword vectors.');
    end
    if length(trblklen) ~= length(sbits)
        error('lte:error','The number of output transport block lengths (size of TRBLKLEN) does not match the number of codeword vectors.');
    end
    
    % If the block lengths were parameterized as a cell array then convert to a vector
    if (iscell(trblklen))
        trblklen = cell2mat(trblklen);
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
        if ~isfield(dch,'NLDPCDecIts')
            dch.NLDPCDecIts = 25;
        end
        [chs(:).NLDPCDecIts] = deal(dch.NLDPCDecIts);
        [chs(:).NBG] = chinfo.NBG;
        if isfield(dch,'TBSLBRM')
            [chs(:).TBSLBRM] = deal(dch.TBSLBRM);
        end
    else
        [chs(:).NTurboDecIts] = deal(dch.NTurboDecIts);
        [chs(:).TxScheme] = deal(dch.TxScheme);
        [chs(:).NIR] = chinfo.NIR;
    end
end