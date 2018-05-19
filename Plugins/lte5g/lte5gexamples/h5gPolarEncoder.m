classdef (StrictDefaults)h5gPolarEncoder < matlab.System & matlab.system.mixin.Propagates ...
    & matlab.system.mixin.CustomIcon
%h5gPolarEncoder Encode binary data using a polar encoder
%   POLARENC = h5gPolarEncoder creates a System object, POLARENC, that
%   encodes binary data using a polar encoder.
%
%   POLARENC = h5gPolarEncoder(Name,Value) creates a polar encoder
%   object, POLARENC, with the specified property Name set to the specified
%   Value. You can specify additional name-value pair arguments in any
%   order as (Name1,Value1, ... , NameN,ValueN).
%
%   POLARENC = h5gPolarEncoder(CODEWORDLENGTH,MESSAGELENGTH,FROZENBITS)
%   creates a polar encoder object, POLARENC, with the CodewordLength
%   property set to CODEWORDLENGTH, MessageLength property set to
%   MESSAGELENGTH, and the FrozenBits property set to FROZENBITS.
%
%   Step method syntax:
%
%   Y = step(POLARENC,X) encodes the binary data, X, using the Polar
%   encoding scheme that you specify by the CodewordLength, MessageLength
%   and FrozenBits properties. It returns the encoded data, Y. Both X and Y
%   are column vectors of data type numeric, or logical where X is of
%   length MessageLength and Y is of length CodewordLength. Depending on
%   the value of the InterleaveInput property, the input X is interleaved
%   prior to encoding.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(obj,x) and y = obj(x) are
%   equivalent.
%
%   h5gPolarEncoder methods:
%
%   step     - Perform polar encoding (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create polar encoder object with same property values
%   isLocked - Locked status (logical)
%
%   h5gPolarEncoder properties:
%
%   CodewordLength    - Code word length (N)
%   MessageLength     - Message length (K)
%   FrozenBits        - Frozen bit vector of length N (F)
%   InterleaveInput   - Interleave input
%
%   % Example: Transmit polar-encoded block of data and decode using
%   %          a Successive Cancellation List decoder.
% 
%   K = 132;            % Message length
%   E = 256;            % Encoded length
%   nMax = 9;           % Maximum value of n, Nmax = 2^n
%   F = h5gPolarConstruct(K,E,nMax);    % {0 info, 1 frozen}
%   N = length(F);      % Mother code length
%   crcLen = 0;         % number of CRC bits
%   nVar = 1.5;         % Noise variance
%   L = 8;              % List length
% 
%   % Object constructions
%   polarEnc  = h5gPolarEncoder(N,K,F,'InterleaveInput',false);
%   bpskMod   = comm.BPSKModulator;
%   chan      = comm.AWGNChannel('NoiseMethod','Variance','Variance',nVar);
%   bpskDemod = comm.BPSKDemodulator('DecisionMethod', ...
%               'Approximate log-likelihood ratio', 'Variance',nVar);
%   polarDec  = h5gPolarDecoder(N,K,F,L,crcLen,'DeinterleaveOutput',false); 
% 
%   % Simulate a frame
%   msg    = randi([0 1],K-crcLen,1);   % Generate a random message
%   enc    = polarEnc(msg);             % Polar encode
%   mod    = bpskMod(enc);              % Modulate
%   rSig   = chan(mod);                 % Add WGN
%   rxLLR  = bpskDemod(rSig);           % Soft demodulate
%   rxBits = polarDec(rxLLR);           % Polar decode
% 
%   % Get bit errors
%   numBitErrs = biterr(rxBits(1:K-crcLen), msg);
%   disp(['Number of bit errors: ' num2str(numBitErrs)]);        
% 
%   See also h5gPolarConstruct, h5gPolarDecoder.

%   Copyright 2017-2018 The MathWorks, Inc.

%   References:
%   [1] 3GPP TS 38.212, "3rd Generation Partnership Project; Technical 
%   Specification Group Radio Access Network; NR; Multiplexing and channel 
%   coding (Release 15), v15.0.0, 2017-12. Sections 5.3.1.

%#codegen

properties (Nontunable)
    %CodewordLength Code word length
    %   Specify the code word length as a positive integer scalar which
    %   also must be a power of 2. The default is 512.
    CodewordLength = 512;
    %MessageLength Message length
    %   Specify the message length as a positive integer scalar. The
    %   default is 256.    
    MessageLength = 256;
    %FrozenBits Frozen bit vector
    %   Specify the frozen bit vector as a binary vector of length
    %   CodewordLength, where a 1-value specifies a frozen bit index and a
    %   0-value specifies the information bit index. The number of 0's in
    %   the vector must be equal to MessageLength. The default is
    %   [ones(256,1);zeros(256,1)].
    %
    %   See also h5gPolarConstruct.
    FrozenBits = [ones(256,1);zeros(256,1)];
end

properties (Nontunable, Logical)
    %InterleaveInput Interleave input 
    %   Specify whether the input bits are interleaved prior to
    %   encoding. The default is true, which corresponds to an iIL value
    %   of 1. The interleaving pattern is as per Section 5.3.1.1.
    InterleaveInput = true;
end

properties(Access = private)
    pG;             % GN matrix
    piInterl;       % interleaver pattern
end

methods
    % CONSTRUCTOR
    function obj = h5gPolarEncoder(varargin)
        setProperties(obj, nargin, varargin{:}, 'CodewordLength',...
            'MessageLength', 'FrozenBits');
    end
    
    function set.CodewordLength(obj, value)
        validateattributes(value, {'numeric'}, ...
            {'real', 'finite', 'positive', 'integer', 'scalar'}, ...
            '', 'CodewordLength'); 

        obj.CodewordLength = value;
    end
    
    function set.MessageLength(obj, value)
        validateattributes(value, {'numeric'}, ...
            {'real', 'finite', 'positive', 'integer', 'scalar'}, ...
            '', 'MessageLength'); 

        obj.MessageLength = value;
    end

    function set.FrozenBits(obj, value)
        validateattributes(value, {'numeric'}, ...
            {'nonempty', 'nonsparse', 'binary', 'vector'}, ...
            '', 'FrozenBits'); 

        obj.FrozenBits = value;
    end

    function set.InterleaveInput(obj,value)
        propName = 'InterleaveInput';
        validateattributes(value, {'logical'}, {'scalar'}, ...
            [class(obj) '.' propName], propName); 

        obj.InterleaveInput = value;
    end
end

methods(Access = protected)   % System object APIs

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
                
        % length(F==0) is K, length(F==1) is N-K
        coder.internal.errorIf(sum(obj.FrozenBits) ~= N-obj.MessageLength, ...
            'lte5g:h5gPolar:InvalidF_NK');
    end

    %% Validate inputs
    function validateInputsImpl(obj, x)
        % Type checks for x
        validateattributes(x, {'double','int8'}, {'nonsparse',...
             'finite', '2d'}, [class(obj) '.step'], 'X');
        
        % Check length of x and MessageLength match
         coder.internal.errorIf(length(x) ~= obj.MessageLength, ...
            'lte5g:h5gPolar:InvalidX_K');
    end

    %% Number of Inputs
    function num = getNumInputsImpl(~)
        num = 1;
    end

    %% Size propagators
    function flag = isInputSizeLockedImpl(~,~)
        flag = true;
    end
        
    function varargout = isOutputFixedSizeImpl(~)
        varargout = {true};    
    end
    
    function varargout = getOutputSizeImpl(obj)
        varargout = {[obj.CodewordLength, 1]}; 
    end
    
    %% Type propagators
    function varargout = getOutputDataTypeImpl(obj)
        varargout = {propagatedInputDataType(obj, 1)};  
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
    
    %% Setup
    function setupImpl(obj,varargin)       
        n = log2(obj.CodewordLength);

        % Get G, nth Kronecker power of kernel 
        ak0 = [1 0; 1 1];   % Arikan's kernel
        allG = cell(n,1);
        allG{1} = ak0;
        for i = 1:n-1
            allG{i+1} = kron(allG{i},ak0);
        end
        obj.pG = allG{n};
        
        % Interleaving pattern
        if obj.InterleaveInput
            obj.piInterl = nr5g.internal.polarInterleaveMap(obj.MessageLength);
        else
            % No interleaving
            obj.piInterl = (0:obj.MessageLength-1)';
        end
    end

    %% Step
    function y = stepImpl(obj,x)

        % Interleave input
        xil = x(obj.piInterl+1);
        
        % Generate u, for only npc==0
        u = zeros(obj.CodewordLength,1);

        % Set input information bits
        u(obj.FrozenBits==0) = xil;

        % Encode using matrix multiplication
        y = mod(u'*obj.pG,2)';  

    end
    
    %%      
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@matlab.System(obj);
        if isLocked(obj)
            s.pG  = obj.pG;
            s.piInterl  = obj.piInterl;
        end
    end
  
    function loadObjectImpl(obj,s,wasLocked)
        if wasLocked
            obj.pG = s.pG;
            obj.piInterl = s.piInterl;
        end
        % Call the base class method
        loadObjectImpl@matlab.System(obj, s);
    end

    function icon = getIconImpl(~)        
        icon = sprintf('Polar\nEncoder');
    end
  
    function varargout = getInputNamesImpl(~)
        % Always have no label for the output port
        varargout = {''};
    end
    
    function varargout = getOutputNamesImpl(~)
        % Always have no label for the output port
        varargout = {''};
    end
    
end

methods(Static, Access=protected)
    
    function header = getHeaderImpl
        header = matlab.system.display.Header('h5gPolarEncoder', ...
            'Title', 'Polar Encoder', 'Text', ...
            'Encode binary data using a polar encoder.');
    end
end

end
