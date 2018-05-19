%h5gRateMatchLDPC 5G LDPC rate matching
%   OUT = h5gRateMatchLDPC(...) rate matches the matrix input data to
%   create the vector OUT. This function includes the stages of bit
%   selection, interleaving defined for LDPC encoded data and code block
%   segmentation (see TS 38.212 Subclauses 5.4.2 and 5.5). The function
%   does not support code block group (re-)transmissions.
%
%   OUT = h5gRateMatchLDPC(IN,OUTLEN,RV,CHS) rate matches the input data IN
%   to create vector OUT of length OUTLEN. The input data is a matrix, each
%   column of which is assumed to be a code block. In the case of the
%   number of columns in IN >= 2, each column is rate matched separately
%   and the results are concatenated into the single output vector OUT. The
%   redundancy version (RV) of the output is controlled by the RV parameter
%   (0,1,2,3). The CHS structure should contain the following fields:
%   NBG        - Base graph number (1,2)
%   Modulation - Modulation format ('BPSK','QPSK','16QAM','64QAM','256QAM')
%   NLayers    - Total number of transmission layers associated with the
%                transport block(s) (1..8)
%   TBSLBRM    - Optional. Transport block size used to configure the soft
%                buffer size for limited buffer rate matching. It is
%                assumed no limit is placed on the number of soft bits,
%                when the field is not provided
%
%   Example:
%   % The example below shows the rate matching of two LDPC encoded code
%   % blocks of length 2000 to a vector of length 5000. The selected 
%   % parameters are: the base graph number is 1, the RV is set to 0, QPSK
%   % modulation is used and the number of layers is 1.
%   
%   % Select channel parameters for rate recovery
%   chs.NBG = 1;
%   chs.rv = 3;
%   chs.Modulation = 'QPSK';
%   chs.NLayers = 1;
% 
%   encoded = ones(6600,2);
%   outlen = 5000;
%   ratematched = h5gRateMatchLDPC(encoded,outlen,chs.rv,chs);
%   size(ratematched)
%
%   See also h5gRateRecoverLDPC, lteRateMatchTurbo.

% Copyright 2017-2018 The MathWorks, Inc.

function out = h5gRateMatchLDPC(in,outlen,rv,chs)
    
    % Ouput empty if the input is empty or the rate matching output length
    % is 0
    if isempty(in) || ~outlen
        out = zeros(0,1,class(in));
        return;
    end

    % Validate the input base graph number
    if ~(chs.NBG == 1 || chs.NBG == 2)
        error('lte:error','The input base graph category number (%d) should be either 1 or 2.',chs.NBG);
    end
    
    % Get modulation order
    if strcmpi(chs.Modulation,'BPSK')
        Qm = 1;
    elseif strcmpi(chs.Modulation,'QPSK')
        Qm = 2;
    elseif strcmpi(chs.Modulation,'16QAM')
        Qm = 4;
    elseif strcmpi(chs.Modulation,'64QAM')
        Qm = 6;
    elseif strcmpi(chs.Modulation,'256QAM')
        Qm = 8;
    end
    
    % Get code block soft buffer size
    N = size(in, 1);
    C = size(in, 2);
    if isfield(chs,'TBSLBRM')
        Nref = floor(3*chs.TBSLBRM/4);
        Ncb = min(N,Nref);
    else % No limit on buffer size
        Ncb = N; 
    end

    % Get starting position in circulr buffer
    if chs.NBG == 1
        Zc = N/66;
        if rv == 0
            k0 = 0;
        elseif rv == 1
            k0 = floor(17*Ncb/(66*Zc))*Zc;
        elseif rv == 2
            k0 = floor(33*Ncb/(66*Zc))*Zc;
        elseif rv == 3
            k0 = floor(56*Ncb/(66*Zc))*Zc;
        end
    else
        Zc = N/50;
        if rv == 0
            k0 = 0;
        elseif rv == 1
            k0 = floor(13*Ncb/(50*Zc))*Zc;
        elseif rv == 2
            k0 = floor(25*Ncb/(50*Zc))*Zc;
        elseif rv == 3
            k0 = floor(43*Ncb/(50*Zc))*Zc;
        end
    end
    
    % Validate the input data size
    if fix(Zc) ~= Zc
        if chs.NBG == 1
            Ncw = 66;
        else
            Ncw = 50;
        end
        error('lte:error','The number of rows in the input data matrix (%d) should be a multiple integer of %d (for base graph %d).', N, Ncw, chs.NBG);
    end
    
    % Get rate matching output for all code blocks and perform code block
    % concatenation according to 38.212 5.4.2 and 5.5
    out = [];
    for r = 0:C-1
        if r <= C-mod(outlen/(chs.NLayers*Qm),C)-1
            E = chs.NLayers*Qm*floor(outlen/(chs.NLayers*Qm*C));
        else
            E = chs.NLayers*Qm*ceil(outlen/(chs.NLayers*Qm*C));
        end
        out = [out; cbsRateMatch(in(:,r+1),E,k0,Ncb,Qm)]; %#okARGOW
    end
    
end

% Rate matching for a single code block segment according to 38.212 5.4.2
function e = cbsRateMatch(d,E,k0,Ncb,Qm)

    % Perform bit selection according to 38.212 5.4.2.1
    k = 0;
    j = 0;
    e = zeros(E,1,class(d));
    while k < E
        if d(mod(k0+j,Ncb)+1) ~= -1 % <NULL> filler bits
            e(k+1) = d(mod(k0+j,Ncb)+1);
            k = k+1;
        end
        j = j+1;
    end
    
    % Perform bit interleaving according to 38.212 5.4.2.2
    e = reshape(e,E/Qm,Qm);
    e = e.';
    e = e(:);
    
end