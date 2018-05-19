classdef (StrictDefaults)h5gPolarDecoder < matlab.System & matlab.system.mixin.Propagates & ...
        matlab.system.mixin.CustomIcon
    %h5gPolarDecoder Decode input using a polar decoder
    %   POLARDEC = h5gPolarDecoder creates a polar decoder System object,
    %   POLARDEC to decode a polar-code encoded bit stream. This object
    %   uses the CRC-Aided Successive-Cancellation List (CA-SCL) decoding
    %   algorithm.
    %
    %   POLARDEC = h5gPolarDecoder(Name,Value) creates a polar
    %   decoder object, POLARDEC, with the specified property Name set to
    %   the specified Value. You can specify additional name-value pair
    %   arguments in any order as (Name1,Value1, ... ,NameN,ValueN).
    %
    %   POLARDEC = h5gPolarDecoder(CODEWORDLENGTH,MESSAGELENGTH, ...
    %   FROZENBITS,LISTLENGTH,CRCLENGTH) creates a polar decoder object,
    %   POLARDEC, with the CodewordLength property set to CODEWORDLENGTH,
    %   MessageLength property set to MESSAGELENGTH, FrozenBits property
    %   set to FROZENBITS, the ListLength property set to LISTLENGTH, and
    %   the CRCLength property set to CRCLENGTH.
    %
    %   Step method syntax:
    %
    %   Y = step(POLARDEC,INPLLR) decodes the received bits, input as
    %   log-likelihood ratios, INPLLR, using the CA-SCL decoding algorithm.
    %   The input is:
    %       INPLLR: a N-by-1 double column containing log-likelihood ratios
    %   for each of the received bits.
    %   The output is:
    %       Y: a double column vector containing the decoded bits. The size
    %       of the output is K, where K represents the number of
    %       information bits encoded.
    %
    %   System objects may be called directly like a function instead of
    %   using the step method. For example, y = step(obj,x) and y = obj(x)
    %   are equivalent.
    %
    %   h5gPolarDecoder methods:
    %
    %   step     - Perform polar decoding (see above)
    %   release  - Allow property value and input characteristics changes
    %   clone    - Create polar decoder object with same property values
    %   isLocked - Locked status (logical)
    %
    %   h5gPolarDecoder properties:
    %
    %   CodewordLength     - Code word length (N)
    %   MessageLength      - Message length (K)
    %   FrozenBits         - Frozen bit vector of length N (F)
    %   ListLength         - List decoding length (L)
    %   CRCLength          - Number of CRC bits
    %   DeinterleaveOutput - Deinterleave output 
    %
    %   % Example: Transmit polar-encoded block of data and decode using
    %   %          a CRC-Aided Successive Cancellation List decoder.
    %
    %   K = 132;            % Message length
    %   E = 256;            % Encoded length
    %   nMax = 9;           % Maximum value of n, Nmax = 2^n
    %   F = h5gPolarConstruct(K,E,nMax);    % {0 info, 1 frozen}
    %   N = length(F);      % Mother code length
    %   crcLen = 24;        % Number of CRC bits
    %   nVar = 1.5;         % Noise variance
    %   L = 8;              % List length
    %
    %   % Object constructions
    %   polarEnc  = h5gPolarEncoder(N,K,F,'InterleaveInput',true);
    %   bpskMod   = comm.BPSKModulator;
    %   chan      = comm.AWGNChannel('NoiseMethod','Variance','Variance',nVar);
    %   bpskDemod = comm.BPSKDemodulator('DecisionMethod', ...
    %               'Approximate log-likelihood ratio', 'Variance',nVar);
    %   polarDec  = h5gPolarDecoder(N,K,F,L,crcLen,'DeinterleaveOutput',true); 
    %
    %   % Simulate a frame
    %   msg    = randi([0 1],K-crcLen,1);   % Generate a random message
    %   msgcrc = h5gCRCEncode(msg,'24C');   % CRC encode
    %   enc    = polarEnc(msgcrc);          % Polar encode
    %   mod    = bpskMod(enc);              % Modulate
    %   rSig   = chan(mod);                 % Add WGN
    %   rxLLR  = bpskDemod(rSig);           % Soft demodulate
    %   rxBits = polarDec(rxLLR);           % Polar decode
    %
    %   % Get bit errors
    %   numBitErrs = biterr(rxBits(1:K-crcLen), msg);
    %   disp(['Number of bit errors: ' num2str(numBitErrs)]);        
    %
    %   See also h5gPolarEncoder, h5gPolarConstruct.
    
    %   Copyright 2017-2018 The MathWorks, Inc.
    
    %#codegen
    
    % References:
    % [1] Tal, I, and Vardy, A., "List decoding of Polar Codes", IEEE
    % Transactions on Information Theory, vol. 61, No. 5, pp. 2213-2226,
    % May 2015. 
    % [2] Stimming, A. B., Parizi, M. B., and Burg, A., "LLR-Based
    % Successive Cancellation List Decoding of Polar Codes", IEEE
    % Transaction on Signal Processing, vol. 63, No. 19, pp.5165-5179,
    % 2015.
    
    properties (Nontunable)
        %CodewordLength Code word length
        %   Specify the code word length as a positive integer scalar which
        %   also must be a power of 2. The default is 512.
        CodewordLength = 512;
        %MessageLength Message length
        %   Specify the message length as a positive integer scalar. The
        %   default is 256.
        MessageLength = 256;
    end

    properties (Nontunable)
        %FrozenBits Frozen bit vector
        %   Specify the frozen bit vector as a binary vector of length
        %   CodewordLength, where a 1-value specifies a frozen bit index
        %   and a 0-value specifies the information bit index. The number
        %   of 0's in the vector must be equal to MessageLength. The
        %   default is [ones(256,1);zeros(256,1)].
        %
        %   See also h5gPolarConstruct.
        FrozenBits = [ones(256,1);zeros(256,1)];
    end

    properties (Nontunable)
        %ListLength List decoding length
        %   Specify the list decoding length as a positive integer scalar.
        %   The default is 4.
        ListLength = 4;     % List decoding length
        %CRCLength Number of CRC bits used for decoding
        %   Specify the number of CRC bits used for decoding. A 0-value
        %   implies no CRC-aided decoding. When enabled, the output decoded
        %   bits include the CRC bits. Supported values include one of [0 6
        %   11 24]. The default is 0.
        CRCLength = 0;      % CRC length
    end
    
    properties (Nontunable, Logical)
        %DeinterleaveOutput Deinterleave output 
        %   Specify whether the decoded bits are deinterleaved prior to
        %   output. The default is true, which corresponds to an iIL value
        %   of 1.
        DeinterleaveOutput = true;
    end

    properties (Access=private)
        % Tal & Vardy data structures
        %   CodewordLength: N, MessageLength: K, FrozenBits: F
        %   ListLength: L
        m;                      % log2(N)
        br;                     % bit reversed vector (N-by-1)
        arrayPtrLLR;            % (m+1)-by-L
        arrayPtrC;              % (m+1)-by-L
        llrPathMetric;          % L-by-1
        pathIdxToArrayIdx;      % (m+1)-by-L
        inactiveArrayIndices;   % (m+1)-by-L
        inactiveArrayIndicesLen;% (m+1)-by-1, (m+1) stack depths
        arrayReferenceCount;    % (m+1)-by-L
        inactivePathIndices;    % L-by-1
        inactivePathIndicesLen; % 1-by-1, stack depth
        activePath;             % L-by-1
        savedCWs;               % N-by-L
        crcDet;                 % CRC-detector
        piInterl;               % interleaver pattern
    end
    
    
    methods
        % Class constructor
        function obj = h5gPolarDecoder(varargin)
            setProperties(obj, nargin, varargin{:}, 'CodewordLength', ...
                'MessageLength', 'FrozenBits', 'ListLength', 'CRCLength');
        end
        
        % Set class parameters
        function set.CodewordLength(obj,value)
            validateattributes(value, {'numeric'}, ...
                {'real', 'finite', 'positive', 'integer', 'scalar'}, ...
                '', 'CodewordLength');
            
            obj.CodewordLength = value;
        end
        
        function set.MessageLength(obj,value)
            validateattributes(value, {'numeric'}, ...
                {'real', 'finite', 'positive', 'integer', 'scalar'}, ...
                '', 'MessageLength');
            
            obj.MessageLength = value;
        end
        
        function set.FrozenBits(obj,value)
            validateattributes(value, {'numeric'}, ...
                {'nonempty', 'nonsparse', 'binary', 'vector'}, ...
                '', 'FrozenBits');
            
            obj.FrozenBits = value;
        end
        
        function set.ListLength(obj,value)
            validateattributes(value, {'numeric'}, ...
                {'real', 'finite', 'positive', 'integer', 'scalar'}, ...
                '', 'ListLength');
            
            obj.ListLength = value;
        end
        
        function set.CRCLength(obj,value)
            validateattributes(value, {'numeric'}, ...
                {'real', 'finite', 'nonnegative', 'integer', 'scalar'}, ...
                '', 'CRCLength');
            
            obj.CRCLength = value;
        end
                
        function set.DeinterleaveOutput(obj,value)
            propName = 'DeinterleaveOutput';
            validateattributes(value, {'logical'}, {'scalar'}, ...
                [class(obj) '.' propName], propName); 

            obj.DeinterleaveOutput = value;
        end
        
    end
    
    methods (Access = protected)
        
        %% Validate properties
        function validatePropertiesImpl(obj)
            
            N = obj.CodewordLength;
            
            % K < N
            coder.internal.errorIf(obj.MessageLength>=N, ...
                'lte5g:h5gPolar:InvalidN_K');
            
            % N is power of 2
            coder.internal.errorIf(log2(N)~=floor(log2(N)), ...
                'lte5g:h5gPolar:Invalid_N');
            
            % F is vector of length N
            coder.internal.errorIf(length(obj.FrozenBits)~=N, ...
                'lte5g:h5gPolar:InvalidF_N');
            
            % length(F==0) is K, length(F==1) is N-K,
            coder.internal.errorIf(sum(obj.FrozenBits) ~= N-obj.MessageLength, ...
                'lte5g:h5gPolar:InvalidF_NK');
            
            % CRCLength: only a few values are supported
            coder.internal.errorIf(~ismember(obj.CRCLength, [0 6 11 24]), ...
                'lte5g:h5gPolar:InvalidCRCLength');
            
        end
        
        %% Validate inputs
        function validateInputsImpl(obj,inpLLR)
            % Type checks for x
            validateattributes(inpLLR, {'double'}, {'nonsparse',...
                'finite', '2d'}, [class(obj) '.step'], 'inpLLR');
            
            % Check length of llrIn and CodewordLength match
            coder.internal.errorIf(length(inpLLR) ~= obj.CodewordLength, ...
                'lte5g:h5gPolar:InvalidIn_N');
            
        end
        
        %% Set number of inputs and outputs
        function numIn = getNumInputsImpl(~)
            numIn = 1; % for llr input
        end
        
        function numOut = getNumOutputsImpl(~)
            numOut = 1;
        end
        
        %% Size propagators
        function flag = isInputSizeLockedImpl(~,~)
            flag = false;
        end
        
        function varargout = isOutputFixedSizeImpl(obj)
            varargout = {propagatedInputFixedSize(obj, 1)};
        end
        
        function varargout = getOutputSizeImpl(obj)
            varargout = {[obj.MessageLength 1]};
        end
        
        %% Type propagators
        function varargout = getOutputDataTypeImpl(~)
            varargout = {'double'};
        end
        
        %% Complexity propagators
        function flag = isInputComplexityLockedImpl(~,~)
            flag = true;
        end
        
        function flag = isOutputComplexityLockedImpl(~,~)
            flag = true;
        end
        
        function varargout = isOutputComplexImpl(~)
            varargout = {false};
        end
        
        %% setupImpl
        function setupImpl(obj,varargin)
            obj.m = log2(obj.CodewordLength);
            
            obj.br = zeros(obj.CodewordLength,1);
            for idxBit = 0:obj.CodewordLength-1
                % 0-based indexing
                obj.br(idxBit+1) = h5gPolarBitReverse(idxBit,obj.m);
            end
            
            % CRC's off 38.212, Section 5.1
            switch obj.CRCLength

                case 6 % binary form, for uplink
                    obj.crcDet = comm.CRCDetector('Polynomial', ...
                        h5gCRCPoly('6'));
                    
                case 11 % binary form, for uplink
                    obj.crcDet = comm.CRCDetector('Polynomial', ...
                        h5gCRCPoly('11'));
                                      
                case 24 % binary form, for downlink
                    obj.crcDet = comm.CRCDetector('Polynomial', ...
                        h5gCRCPoly('24C'));
                    
                otherwise % default for codegen, no CRC
                    obj.crcDet = comm.CRCDetector('Polynomial', 1);
            end
            
            % Interleaving pattern
            if obj.DeinterleaveOutput
                obj.piInterl = nr5g.internal.polarInterleaveMap(obj.MessageLength);
            else
                % No deinterleaving
                obj.piInterl = (0:obj.MessageLength-1).';
            end
        end
        
        %% stepImpl
        function dec = stepImpl(obj,llrIn)
            
            % Initialization
            h5gPolarDecoder.initializeDataStructures(obj);
            iniPathIdx = h5gPolarDecoder.assignInitialPath(obj);
            sp = h5gPolarDecoder.getArrayPtrP(obj, 1, iniPathIdx);
            obj.arrayPtrLLR{1,sp}(:,1) = llrIn(obj.br+1);  % LLRs
            mplus1 = obj.m+1;
            
            % Main loop
            for phase = 1:obj.CodewordLength
                if isempty(coder.target)
                    h5gPolarDecoder.recursivelyCalcP(obj, mplus1, phase);
                else
                    % Recursive fcn workaround for codegen
                    sttStr = h5gPolarDecoder.obj2SttStr(obj);
                    sttStr = h5gPolarDecoder.recCalcP(sttStr, mplus1, phase);
                    h5gPolarDecoder.sttStr2obj(obj, sttStr);
                end
                
                pm2 = mod(phase-1,2);
                if obj.FrozenBits(phase)==1
                    % Path for frozen (and punctured) bits
                    for pathIdx = 1:obj.ListLength
                        if ~obj.activePath(pathIdx)
                            continue;
                        end
                        sc = h5gPolarDecoder.getArrayPtrC(obj, mplus1, pathIdx);
                        obj.arrayPtrC{mplus1,sc}(1,pm2+1) = 0; % set to 0
                        
                        obj.llrPathMetric(pathIdx) = ...
                            obj.llrPathMetric(pathIdx) + ...
                            log(1 + exp(-obj.arrayPtrLLR{mplus1,sc}(1,1)));
                    end
                else % Path for info bits
                    h5gPolarDecoder.contPathsUnfrozenBit(obj, phase);
                end
                
                if pm2==1
                    if isempty(coder.target)
                        h5gPolarDecoder.recursivelyUpdateC(obj, mplus1, phase);
                    else
                        % Recursive fcn workaround for codegen
                        sttStr = h5gPolarDecoder.obj2SttStr(obj);
                        sttStr = h5gPolarDecoder.recUpdateC(sttStr, mplus1, phase);
                        h5gPolarDecoder.sttStr2obj(obj, sttStr);
                    end
                end
            end
            
            % Return the best codeword in the list
            %   Use CRC checks, if enabled
            pathIdx1 = 1;
            p1 = realmax;
            crcCW = false;
            for pathIdx = 1:obj.ListLength
                if ~obj.activePath(pathIdx)
                    continue;
                end
                
                sc = h5gPolarDecoder.getArrayPtrC(obj, mplus1, pathIdx);
                if obj.CRCLength>0
                    canCW = obj.savedCWs(:,sc);     % N, with frozen bits
                    canMsg = canCW(obj.FrozenBits==0,1);    % K bits only
                    canMsg(obj.piInterl+1) = canMsg;    % deinterleave
                    if h5gPolarDecoder.crcCheck(obj,canMsg) % 1 => fail
                        continue; % move to next path
                    end
                end
                crcCW = true;
                if p1 > obj.llrPathMetric(pathIdx)
                    p1 = obj.llrPathMetric(pathIdx);
                    pathIdx1 = pathIdx;
                end
            end
            
            if ~crcCW   % no codeword found which passes crcCheck
                pathIdx1 = 1;
                p1 = realmax;
                for pathIdx = 1:obj.ListLength
                    if ~obj.activePath(pathIdx)
                        continue;
                    end
                    
                    if p1 > obj.llrPathMetric(pathIdx)
                        p1 = obj.llrPathMetric(pathIdx);
                        pathIdx1 = pathIdx;
                    end
                end
            end
            
            % Get decoded bits
            sc = h5gPolarDecoder.getArrayPtrC(obj, mplus1, pathIdx1);
            decCW = obj.savedCWs(:,sc);         % N, with frozen bits
            dec = decCW(obj.FrozenBits==0,1);   % K, info bits only
            
            % Deinterleave output
            dec(obj.piInterl+1) = dec;
            
        end
               
        %%
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@matlab.System(obj);
            if isLocked(obj)
                s.m      = obj.m;
                s.br     = obj.br;
                s.crcDet = obj.crcDet;
                s.piInterl = obj.piInterl;
                % Others not added as they get initialized
            end
        end
        
        function loadObjectImpl(obj, s, wasLocked)
            if wasLocked
                obj.m      = s.m;
                obj.br     = s.br;
                obj.crcDet = s.crcDet;
                obj.piInterl = s.piInterl;
                % Others not added as they get initialized
            end
            % Call the base class method
            loadObjectImpl@matlab.System(obj, s);
        end
        
        function icon = getIconImpl(~)
            icon = sprintf('Polar\nDecoder');
        end
        
        function varargout = getInputNamesImpl(~)
            % Always have no label for the input port (only one)
            varargout = {''};
        end
        
        function varargout = getOutputNamesImpl(~)
            % Always have no label for the output port (only one)
            varargout = {''};
        end
        
    end
    
    methods (Static, Access = protected)
        function header = getHeaderImpl
            header = matlab.system.display.Header('h5gPolarDecoder', ...
                'Title', 'Polar Decoder');
        end
    end
    
    methods (Static, Access = private)
        
        function errFlag = crcCheck(obj, canMsg)
            if obj.CRCLength>0
                % errFlag is 1 for error, 0 for no errors
                [~, errFlag] = obj.crcDet(canMsg);
            end
        end
        
        %-----------------------------------------------------------------
        function contPathsUnfrozenBit(obj, phase)
            % Input:
            %   phase: phase phi, 1-based, 1:2^m, or 1:N
            
            % Populate probForks
            probForks = -realmax*ones(obj.ListLength,2);
            i = 0;
            mplus1 = obj.m+1;
            for pathIdx = 1:obj.ListLength
                if obj.activePath(pathIdx)
                    sp = h5gPolarDecoder.getArrayPtrP(obj, mplus1, pathIdx);
                    probForks(pathIdx,1) = - (obj.llrPathMetric(pathIdx) ...
                        + log(1+exp(-(obj.arrayPtrLLR{mplus1,sp}(1,1)))));
                    probForks(pathIdx,2) = - (obj.llrPathMetric(pathIdx) ...
                        + log(1+exp((obj.arrayPtrLLR{mplus1,sp}(1,1)))));
                    i = i+1;
                end
            end
            
            rho = min(2*i,obj.ListLength);
            contForks = zeros(obj.ListLength,2);
            % Populate contForks such that contForks(l,b) is true iff
            % probForks(l,b) is one of rho largest entries in probForks.
            prob = sort(probForks(:), 'descend');
            if rho>0
                threshold = prob(rho);
            else
                threshold = prob(1); % Largest
            end
            numPop = 0;
            for pathIdx = 1:obj.ListLength
                for bIdx = 1:2
                    if numPop==rho
                        break;
                    end
                    if probForks(pathIdx,bIdx)>threshold
                        contForks(pathIdx,bIdx) = 1;
                        numPop = numPop+1;
                    end
                end
            end
            
            if numPop<rho
                for pathIdx = 1:obj.ListLength
                    for bIdx = 1:2
                        if numPop==rho
                            break;
                        end
                        if probForks(pathIdx,bIdx)==threshold
                            contForks(pathIdx,bIdx) = 1;
                            numPop = numPop+1;
                        end
                    end
                end
            end
            
            % First, kill-off non-continuing paths
            for pathIdx = 1:obj.ListLength
                if ~obj.activePath(pathIdx)
                    continue;
                end
                if contForks(pathIdx,1)==0 && contForks(pathIdx,2)==0
                    h5gPolarDecoder.killPath(obj, pathIdx);
                end
            end
            
            % Continue relevant paths, duplicating if necessary
            pm2 = mod(phase-1,2);
            for pathIdx = 1:obj.ListLength
                if contForks(pathIdx,1)==0 && contForks(pathIdx,2)==0
                    % Both forks are bad
                    continue;
                end
                
                sc = h5gPolarDecoder.getArrayPtrC(obj, mplus1, pathIdx);
                if contForks(pathIdx,1)==1 && contForks(pathIdx,2)==1
                    % Both forks are good
                    obj.arrayPtrC{mplus1,sc}(1,pm2+1) = 0;
                    obj.savedCWs(phase,sc) = 0;
                    
                    pathIdx1 = h5gPolarDecoder.clonePath(obj, pathIdx);
                    sc2 = h5gPolarDecoder.getArrayPtrC(obj, mplus1, pathIdx1);
                    obj.savedCWs(1:phase-1,sc2) = obj.savedCWs(1:phase-1,sc);
                    
                    obj.arrayPtrC{mplus1,sc2}(1,pm2+1) = 1;
                    obj.savedCWs(phase,sc2) = 1;
                    obj.llrPathMetric(pathIdx) = obj.llrPathMetric(pathIdx) ...
                        + log(1 + exp(-obj.arrayPtrLLR{mplus1,sc}(1,1)));
                    obj.llrPathMetric(pathIdx1) = obj.llrPathMetric(pathIdx1) ...
                        + log(1 + exp(obj.arrayPtrLLR{mplus1,sc2}(1,1)));
                else
                    % Exactly one fork is good
                    if contForks(pathIdx,1)==1
                        obj.arrayPtrC{mplus1,sc}(1,pm2+1) = 0;
                        obj.savedCWs(phase,sc) = 0;
                        obj.llrPathMetric(pathIdx) = obj.llrPathMetric(pathIdx) ...
                            + log(1 + exp(-obj.arrayPtrLLR{mplus1,sc}(1,1)));
                    else
                        obj.arrayPtrC{mplus1,sc}(1,pm2+1) = 1;
                        obj.savedCWs(phase,sc) = 1;
                        obj.llrPathMetric(pathIdx) = obj.llrPathMetric(pathIdx) ...
                            + log(1 + exp(obj.arrayPtrLLR{mplus1,sc}(1,1)));
                    end
                end
            end
            
        end
        
        %----------------------------------------------------------------
        % Mid-level fcns
        %----------------------------------------------------------------
        function recursivelyCalcP(obj, layer, phase)
            % Input:
            %   layer: layer lambda, 1-based, 1:m+1
            %   phase: phase phi, 1-based, 1:2^layer or 1:N
            % Only for simulation
            
            if layer==1
                return;
            end
            psi = floor((phase-1)/2)+1;
            pm2 = mod(phase-1,2);
            if pm2==0
                h5gPolarDecoder.recursivelyCalcP(obj, layer-1, psi);
            end
            
            expm = 2^(obj.m-layer+1);
            for pathIdx = 1:obj.ListLength
                if ~obj.activePath(pathIdx)
                    continue;
                end
                sp = h5gPolarDecoder.getArrayPtrP(obj, layer, pathIdx);
                spminus1 = h5gPolarDecoder.getArrayPtrP(obj, layer-1, pathIdx);
                sc = h5gPolarDecoder.getArrayPtrC(obj, layer, pathIdx);
                
                for beta = 0:expm-1
                    if pm2==0
                        % LLR
                        aa = obj.arrayPtrLLR{layer-1,spminus1}( 2*beta+1,1 );
                        bb = obj.arrayPtrLLR{layer-1,spminus1}( 2*beta+2,1 );
                        % Consider a table look-up instead
                        if max( abs(aa), abs(bb)) < 40
                            % Equation 8a, Stimming
                            obj.arrayPtrLLR{layer,sp}(beta+1,1) = ...
                                log( (exp(aa+bb)+1)/(exp(aa) + exp(bb)) );
                        else
                            % Equation 9, Stimming
                            obj.arrayPtrLLR{layer,sp}(beta+1,1) = ...
                                sign(aa)*sign(bb)*min(abs(aa), abs(bb));
                        end
                    else
                        u1 = obj.arrayPtrC{layer,sc}(beta+1,1);
                        % Equation 8b, Stimming
                        obj.arrayPtrLLR{layer,sp}(beta+1,1) = (-1)^u1 * ...
                            obj.arrayPtrLLR{layer-1,spminus1}( 2*beta+1,1 ) + ...
                            obj.arrayPtrLLR{layer-1,spminus1}( 2*beta+2,1 );
                    end
                end
            end
        end
        
        %-----------------------------------------------------------------
        function sttStr = recCalcP(sttStr, layer, phase)
            % Non-object form of recursivelyCalcP method, still recursive
            % Input:
            %   layer: layer lambda, 1-based, 1:m+1
            %   phase: phase phi, 1-based, 1:2^layer or 1:N
            % Only for codegen
            
            if layer==1
                return;
            end
            psi = floor((phase-1)/2)+1;
            pm2 = mod(phase-1,2);
            if pm2==0
                sttStr = h5gPolarDecoder.recCalcP(sttStr, layer-1, psi);
            end
            
            expm = 2^(sttStr.m-layer+1);
            for pathIdx = 1:sttStr.L
                if ~sttStr.activePath(pathIdx)
                    continue;
                end
                
                [sp, sttStr] = h5gPolarDecoder.getStrArrayPtrP(sttStr, layer, pathIdx);
                [spminus1, sttStr] = h5gPolarDecoder.getStrArrayPtrP(sttStr, layer-1, pathIdx);
                [sc, sttStr] = h5gPolarDecoder.getStrArrayPtrC(sttStr, layer, pathIdx);
                for beta = 0:expm-1
                    if pm2==0
                        % LLR
                        aa = sttStr.arrayPtrLLR{layer-1,spminus1}( 2*beta+1,1 );
                        bb = sttStr.arrayPtrLLR{layer-1,spminus1}( 2*beta+2,1 );
                        % Consider a table look-up instead
                        if max( abs(aa), abs(bb)) < 40
                            % Equation 8a, Stimming
                            sttStr.arrayPtrLLR{layer,sp}(beta+1,1) = ...
                                log( (exp(aa+bb)+1)/(exp(aa) + exp(bb)) );
                        else
                            % Equation 9, Stimming
                            sttStr.arrayPtrLLR{layer,sp}(beta+1,1) = ...
                                sign(aa)*sign(bb)*min(abs(aa), abs(bb));
                        end
                    else
                        u1 = sttStr.arrayPtrC{layer,sc}(beta+1,1);
                        % Equation 8b, Stimming
                        sttStr.arrayPtrLLR{layer,sp}(beta+1,1) = (-1)^u1 * ...
                            sttStr.arrayPtrLLR{layer-1,spminus1}( 2*beta+1,1 ) + ...
                            sttStr.arrayPtrLLR{layer-1,spminus1}( 2*beta+2,1 );
                    end
                end
            end
        end
        
        %-----------------------------------------------------------------
        function recursivelyUpdateC(obj, layer, phase)
            % Input:
            %   layer: layer lambda, 1-based, 1:m+1
            %   phase: phase phi, 1-based, 1:2^layer or 1:N, must be odd
            % Only for simulation
            
            psi = floor((phase-1)/2);
            pm2 = mod(psi,2);
            expm = 2^(obj.m-layer+1);
            for pathIdx = 1:obj.ListLength
                if ~obj.activePath(pathIdx)
                    continue;
                end
                sc = h5gPolarDecoder.getArrayPtrC(obj, layer, pathIdx);
                scminus1 = h5gPolarDecoder.getArrayPtrC(obj, layer-1, pathIdx);
                for beta = 0:expm-1
                    obj.arrayPtrC{layer-1,scminus1}(2*beta+1,pm2+1) = ...
                        xor(obj.arrayPtrC{layer,sc}(beta+1,1), ...
                        obj.arrayPtrC{layer,sc}(beta+1,2));
                    obj.arrayPtrC{layer-1,scminus1}(2*beta+2,pm2+1) = ...
                        obj.arrayPtrC{layer,sc}(beta+1,2);
                end
            end
            
            if pm2==1
                h5gPolarDecoder.recursivelyUpdateC(obj, layer-1, psi+1);
            end
            
        end
        
        %-----------------------------------------------------------------
        function sttStr = recUpdateC(sttStr, layer, phase)
            % Non-object form of recursivelyUpdateC method, still recursive
            % Input:
            %   layer: layer lambda, 1-based, 1:m+1
            %   phase: phase phi, 1-based, 1:2^layer or 1:N, must be odd
            
            psi = floor((phase-1)/2);
            pm2 = mod(psi,2);
            expm = 2^(sttStr.m-layer+1);
            for pathIdx = 1:sttStr.L
                if ~sttStr.activePath(pathIdx)
                    continue;
                end
                [sc, sttStr] = h5gPolarDecoder.getStrArrayPtrC(sttStr, layer, pathIdx);
                [scminus1, sttStr] = h5gPolarDecoder.getStrArrayPtrC(sttStr, layer-1, pathIdx);
                for beta = 0:expm-1
                    sttStr.arrayPtrC{layer-1,scminus1}(2*beta+1,pm2+1) = ...
                        xor(sttStr.arrayPtrC{layer,sc}(beta+1,1), ...
                        sttStr.arrayPtrC{layer,sc}(beta+1,2));
                    sttStr.arrayPtrC{layer-1,scminus1}(2*beta+2,pm2+1) = ...
                        sttStr.arrayPtrC{layer,sc}(beta+1,2);
                end
            end
            
            if pm2==1
                sttStr = h5gPolarDecoder.recUpdateC(sttStr, layer-1, psi+1);
            end
            
        end
        
        %-------------------------------------------------------------------------
        % Low-level fcns
        %-------------------------------------------------------------------------
        function initializeDataStructures(obj)
            % Indices are all 1-based MATLAB indices.
            
            mplus1 = obj.m+1;
            parrayPtrLLR = cell(mplus1,obj.ListLength); % store arrays
            parrayPtrC = cell(mplus1,obj.ListLength); % store arrays
            for layer = 1:mplus1
                expm = 2^(mplus1-layer);
                for s = 1:obj.ListLength
                    parrayPtrLLR{layer,s} = zeros(expm,1);
                    parrayPtrC{layer,s} = zeros(expm,2); % binary-valued: 0,1
                end
            end
            % An extra layer of in-direction helps codegen, update post-17B
            obj.arrayPtrLLR = parrayPtrLLR;
            obj.llrPathMetric = zeros(obj.ListLength,1);
            
            obj.arrayPtrC = parrayPtrC;
            
            obj.pathIdxToArrayIdx = ones(mplus1,obj.ListLength);
            
            obj.inactiveArrayIndices = zeros(mplus1,obj.ListLength);
            obj.inactiveArrayIndicesLen = zeros(mplus1,1);
            for layer = 1:mplus1
                obj.inactiveArrayIndices(layer,:) = 1:obj.ListLength;
                obj.inactiveArrayIndicesLen(layer,1) = obj.ListLength;
            end
            obj.arrayReferenceCount = zeros(mplus1,obj.ListLength);
            
            obj.inactivePathIndices = (1:obj.ListLength).';     % all paths are inactive
            obj.inactivePathIndicesLen = obj.ListLength;
            obj.activePath = zeros(obj.ListLength,1,'logical'); % no paths are active
            
            obj.savedCWs = zeros(obj.CodewordLength,obj.ListLength); % saved codewords
        end
        
        %----------------------------------------------------------------
        function pathIdx = assignInitialPath(obj)
            % Output:
            %   pathIdx: initial path index l, 1-based, 1:L
            
            pathIdx = obj.inactivePathIndices(obj.inactivePathIndicesLen,1);
            obj.inactivePathIndicesLen = obj.inactivePathIndicesLen-1;
            obj.activePath(pathIdx) = true;
            
            % Associate arrays with path index
            for layer = 1:obj.m+1
                s = obj.inactiveArrayIndices(layer,obj.inactiveArrayIndicesLen(layer));
                obj.inactiveArrayIndicesLen(layer) = obj.inactiveArrayIndicesLen(layer)-1;
                
                obj.pathIdxToArrayIdx(layer,pathIdx) = s;
                obj.arrayReferenceCount(layer,pathIdx) = 1;
            end
            
        end
        
        %----------------------------------------------------------------
        function clPathIdx = clonePath(obj, pathIdx)
            % Input:
            %   pathIdx: path index to clone, l, 1-based, 1:L
            % Output:
            %   clPathIdx: cloned path index, l', 1-based
            
            clPathIdx = obj.inactivePathIndices(obj.inactivePathIndicesLen,1);
            obj.inactivePathIndicesLen = obj.inactivePathIndicesLen-1;
            obj.activePath(clPathIdx) = true;
            obj.llrPathMetric(clPathIdx) = obj.llrPathMetric(pathIdx);
            
            % Make clPathIdx (l') reference same arrays as pathIdx (l)
            for layer = 1:obj.m+1
                s = obj.pathIdxToArrayIdx(layer,pathIdx);
                
                obj.pathIdxToArrayIdx(layer,clPathIdx) = s;
                obj.arrayReferenceCount(layer,s) = ...
                    obj.arrayReferenceCount(layer,s)+1;
            end
            
        end
        
        %----------------------------------------------------------------
        function killPath(obj, pathIdx)
            % Input:
            %   pathIdx: path index to kill, l, 1-based, 1:L
            
            % Mark path pathIdx as inactive
            obj.activePath(pathIdx) = false;
            obj.inactivePathIndicesLen = obj.inactivePathIndicesLen+1;
            obj.inactivePathIndices(obj.inactivePathIndicesLen,1) = pathIdx;
            obj.llrPathMetric(pathIdx) = 0;
            
            % Disassociate arrays with path Idx
            for layer = 1:obj.m+1
                s = obj.pathIdxToArrayIdx(layer,pathIdx);
                obj.arrayReferenceCount(layer,s) = ...
                    obj.arrayReferenceCount(layer,s)-1;
                
                if obj.arrayReferenceCount(layer,s)==0
                    obj.inactiveArrayIndicesLen(layer,1) = obj.inactiveArrayIndicesLen(layer,1)+1;
                    obj.inactiveArrayIndices(layer, obj.inactiveArrayIndicesLen(layer,1)) = s;
                end
            end
            
        end
        
        %----------------------------------------------------------------
        function s2 = getArrayPtrP(obj, layer, pathIdx)
            % Input:
            %   layer:   layer lambda, 1-based, 1:m+1
            %   pathIdx: path index l, 1-based, 1:L
            % Output:
            %   s2: corresponding pathIdx for same layer
            
            s = obj.pathIdxToArrayIdx(layer,pathIdx);
            if obj.arrayReferenceCount(layer,s)==1
                s2 = s;
            else
                s2 = obj.inactiveArrayIndices(layer, obj.inactiveArrayIndicesLen(layer,1));
                obj.inactiveArrayIndicesLen(layer,1) = obj.inactiveArrayIndicesLen(layer,1)-1;
                
                % deep copy
                obj.arrayPtrLLR{layer,s2} = obj.arrayPtrLLR{layer,s};
                obj.arrayPtrC{layer,s2} = obj.arrayPtrC{layer,s};
                
                obj.arrayReferenceCount(layer,s) = ...
                    obj.arrayReferenceCount(layer,s)-1;
                obj.arrayReferenceCount(layer,s2) = 1;
                obj.pathIdxToArrayIdx(layer,pathIdx) = s2;
            end
            
        end
        
        %----------------------------------------------------------------
        function [s2, sttStr] = getStrArrayPtrP(sttStr, layer, pathIdx)
            % Non object form of getArrayPtrP
            % Input:
            %   layer:   layer lambda, 1-based, 1:m+1
            %   pathIdx: path index l, 1-based, 1:L
            % Output:
            %   s2: corresponding pathIdx for same layer
            
            s = sttStr.pathIdxToArrayIdx(layer,pathIdx);
            if sttStr.arrayReferenceCount(layer,s)==1
                s2 = s;
            else
                s2 = sttStr.inactiveArrayIndices(layer, sttStr.inactiveArrayIndicesLen(layer,1));
                sttStr.inactiveArrayIndicesLen(layer,1) = sttStr.inactiveArrayIndicesLen(layer,1)-1;
                
                % deep copy
                sttStr.arrayPtrLLR{layer,s2} = sttStr.arrayPtrLLR{layer,s};
                sttStr.arrayPtrC{layer,s2} = sttStr.arrayPtrC{layer,s};
                
                sttStr.arrayReferenceCount(layer,s) = ...
                    sttStr.arrayReferenceCount(layer,s)-1;
                sttStr.arrayReferenceCount(layer,s2) = 1;
                sttStr.pathIdxToArrayIdx(layer,pathIdx) = s2;
            end
            
        end
        
        %----------------------------------------------------------------
        function s2 = getArrayPtrC(obj, layer, pathIdx)
            % Input:
            %   layer:   layer lambda, 1-based, 1:m+1
            %   pathIdx: path index l, 1-based, 1:L
            % Output:
            %   ptrC: corresponding pathIdx for same layer
            
            s = obj.pathIdxToArrayIdx(layer,pathIdx);
            if obj.arrayReferenceCount(layer,s)==1
                s2 = s;
            else
                s2 = obj.inactiveArrayIndices(layer,obj.inactiveArrayIndicesLen(layer,1));
                obj.inactiveArrayIndicesLen(layer,1) = obj.inactiveArrayIndicesLen(layer,1)-1;
                
                % deep copy
                obj.arrayPtrC{layer,s2} = obj.arrayPtrC{layer,s};
                obj.arrayPtrLLR{layer,s2} = obj.arrayPtrLLR{layer,s};
                obj.arrayReferenceCount(layer,s) = ...
                    obj.arrayReferenceCount(layer,s)-1;
                obj.arrayReferenceCount(layer,s2) = 1;
                obj.pathIdxToArrayIdx(layer,pathIdx) = s2;
            end
            
        end
        
        %----------------------------------------------------------------
        function [s2, sttStr] = getStrArrayPtrC(sttStr, layer, pathIdx)
            % Non object form of getArrayPtrP
            % Input:
            %   layer:   layer lambda, 1-based, 1:m+1
            %   pathIdx: path index l, 1-based, 1:L
            % Output:
            %   ptrC: corresponding pathIdx for same layer
            
            s = sttStr.pathIdxToArrayIdx(layer,pathIdx);
            if sttStr.arrayReferenceCount(layer,s)==1
                s2 = s;
            else
                s2 = sttStr.inactiveArrayIndices(layer,sttStr.inactiveArrayIndicesLen(layer,1));
                sttStr.inactiveArrayIndicesLen(layer,1) = sttStr.inactiveArrayIndicesLen(layer,1)-1;
                
                % deep copy
                sttStr.arrayPtrC{layer,s2} = sttStr.arrayPtrC{layer,s};
                sttStr.arrayPtrLLR{layer,s2} = sttStr.arrayPtrLLR{layer,s};
                
                sttStr.arrayReferenceCount(layer,s) = ...
                    sttStr.arrayReferenceCount(layer,s)-1;
                sttStr.arrayReferenceCount(layer,s2) = 1;
                sttStr.pathIdxToArrayIdx(layer,pathIdx) = s2;
            end
            
        end
        
        %----------------------------------------------------------------
        function sttStr = obj2SttStr(obj)
            % Pass out obj states (private properties) as a structure, sttStr.
            
            % Pack private properties into a state struct, sttStr
            sttStr.L                        = obj.ListLength;
            sttStr.m                        = obj.m;                       % log2(N)
            sttStr.arrayPtrLLR              = obj.arrayPtrLLR;             % (m+1)-by-L
            sttStr.llrPathMetric            = obj.llrPathMetric;           % L-by-1
            sttStr.arrayPtrC                = obj.arrayPtrC;               % (m+1)-by-L
            sttStr.pathIdxToArrayIdx        = obj.pathIdxToArrayIdx;       % (m+1)-by-L
            sttStr.inactiveArrayIndices     = obj.inactiveArrayIndices;    % (m+1)-by-L
            sttStr.inactiveArrayIndicesLen  = obj.inactiveArrayIndicesLen; % (m+1)-by-1, (m+1) stack depths
            sttStr.arrayReferenceCount      = obj.arrayReferenceCount;     % (m+1)-by-L
            sttStr.inactivePathIndices      = obj.inactivePathIndices;     % L-by-1
            sttStr.inactivePathIndicesLen   = obj.inactivePathIndicesLen;  % 1-by-1, stack depth
            sttStr.activePath               = obj.activePath;              % L-by-1
            sttStr.savedCWs                 = obj.savedCWs;                % L-by-N
        end
        
        %----------------------------------------------------------------
        function sttStr2obj(obj, sttStr)
            % Update obj states (private properties) from the sttStr
            
            % Update obj with sttStr fields, no need for m, L
            obj.arrayPtrLLR              = sttStr.arrayPtrLLR;
            obj.llrPathMetric            = sttStr.llrPathMetric;
            obj.arrayPtrC                = sttStr.arrayPtrC;
            obj.pathIdxToArrayIdx        = sttStr.pathIdxToArrayIdx;
            obj.inactiveArrayIndices     = sttStr.inactiveArrayIndices;
            obj.inactiveArrayIndicesLen  = sttStr.inactiveArrayIndicesLen;
            obj.arrayReferenceCount      = sttStr.arrayReferenceCount;
            obj.inactivePathIndices      = sttStr.inactivePathIndices;
            obj.inactivePathIndicesLen   = sttStr.inactivePathIndicesLen;
            obj.activePath               = sttStr.activePath;
            obj.savedCWs                 = sttStr.savedCWs;
            
        end
        
    end
    
end

%-------------------------------------------------------------------------
function poly = h5gCRCPoly(tag)
%
%   poly = h5gCRCPoly(TAG) returns the CRC polynomial specified as per 
%   Section 5.1, TS 38.212 for the given TAG. TAG must be one of 
%   {'24C','11','6'} strings.
%
%   % Example: Generate the 24C polynomial
%   p24c = h5gCRCPoly('24C');
% 
%   See also h5gCRCEncode.

%   Reference:
%   [1] 3GPP TS 38.212, "3rd Generation Partnership Project; Technical 
%   Specification Group Radio Access Network; NR; Multiplexing and channel 
%   coding (Release 15), v15.0.0, 2017-12. Section 5.1.

% Use poly as polynomial specification in comm.CRCGenerator
% poly is a binary vector, in descending order
switch tag
    case '24C'  % Polar, downlink
        poly = [1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1];
    case '11'   % Polar, uplink, K>30
        poly = [1 1 1 0 0 0 1 0 0 0 0 1];
    case '6'    % Polar, uplink, 18<=K<=25
        poly = [1 1 0 0 0 0 1];
end

end

% [EOF]
