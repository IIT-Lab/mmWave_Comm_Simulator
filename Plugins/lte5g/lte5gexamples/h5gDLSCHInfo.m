%h5gDLSCHInfo 5G DL-SCH information
%   INFO = h5gDLSCHInfo(CHS,OUTLEN,TRBLKLEN) returns the structure
%   INFO containing the Downlink Shared Channel (DL-SCH) related
%   information.
%
%   OUTLEN is an input vector (one or two elements in length) defining the
%   codeword length(s) that the input transport block(s) should be rate
%   matched to. TRBLKLEN is an input vector (one or two elements in length)
%   defining the transport block sizes.
%
%   INFO contains the following fields:
%   C       - Total number of code blocks
%   NLayers - Number of layers associated with transport block/codeword
%   Qm      - Bits per symbol variable used in rate matching calculation
%   RV      - RV value associated with one codeword (if RV present 
%             at input)
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
%   CHS is an input parameter structure defining aspects of the Physical
%   Downlink Shared Channel (PDSCH) onto which the codeword(s) will be
%   mapped, and the soft buffer size and redundancy version(s) of the
%   generated codeword(s). The required fields are:
%   Modulation     - A string array or cell array specifying the modulation
%                    format for one or two transport blocks
%                    ('QPSK','16QAM','64QAM','256QAM')
%   NLayers        - Total number of transmission layers associated with 
%                    the transport block(s) (1..8)
%   RV             - Vector of 1 or 2 redundancy version indicators
%                    (0,1,2,3)
%   TargetCodeRate - Target coding rate used together with transport block
%                    size to determine the LDPC base graph
%   TBSLBRM        - Optional. Transport block size used to configure the 
%                    soft buffer size for limited buffer rate matching. It 
%                    is assumed no limit is placed on the number of soft 
%                    bits, when the field is not provided

%   Copyright 2018 The MathWorks, Inc.

function info = h5gDLSCHInfo(chs,Gvec,Avec)

    ncw = length(Avec);
    info = struct();
    for icw = 1:ncw
        
        G = Gvec(icw);
        A = Avec(icw);
        
        % Get modulation order
        if iscell(chs.Modulation)
            Modulation = chs.Modulation(mod(icw-1,length(chs.Modulation))+1);
        else
            Modulation = chs.Modulation;
        end
        if strcmpi(Modulation,'QPSK')
            Qm = 2;
        elseif strcmpi(Modulation,'16QAM')
            Qm = 4;
        elseif strcmpi(Modulation,'64QAM')
            Qm = 6;
        elseif strcmpi(Modulation,'256QAM')
            Qm = 8;
        end
        
        % LDPC base graph selection
        if A <= 292 || (A<=3824 && chs.TargetCodeRate<=0.67) || chs.TargetCodeRate<=0.25
            nbg = 2;
        else
            nbg = 1;
        end

        % Get transport block size after CRC attachement according to 38.212
        % 6.2.1 and 7.2.1
        if A > 3824
            L = 24;
        else
            L = 16;
        end
        B = A + L;

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
            C = ceil(B/(Kcb-24));
            Bd = B+C*24;
        end

        % Obtain and verify the number of bits per code block (excluding CB-CRC
        % bits)
        cbz = B/C;
        if fix(cbz) ~= (cbz)
            error('lte:error','The length of the input data (%d) is not a multiple integer of the number of code blocks (%d)',B,C);
        end

        % Number of bits in each code block (excluding Filler bits)
        Kd = Bd/C;

        % Find the minimum value of Z in all sets of lifting sizes in Table
        % 5.3.2-1, denoted as Zc, such that Kb*Zc>=Kd
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
        Zlist = [2:16 18:2:32 36:4:64 72:8:112 114 120 128 160:16:256 288 320 352 384]; 
        for Zc = Zlist
            if Kb*Zc >= Kd
                break;
            end
        end

        % Get number of bits (including Filler bits) to be encoded by LDPC
        % encoder, and length of circular buffer 
        if nbg == 1
            K = 22*Zc;
            N = 66*Zc;
        else
            K = 10*Zc;
            N = 50*Zc;
        end

        % Get rate matching output size for each code blocks
        E = zeros(C,1);
        for r = 0:C-1
            if r <= C-mod(G/(chs.NLayers*Qm),C)-1
                E(r+1) = chs.NLayers*Qm*floor(G/(chs.NLayers*Qm*C));
            else
                E(r+1) = chs.NLayers*Qm*ceil(G/(chs.NLayers*Qm*C));
            end
        end

        % Get code block soft buffer size
        if isfield(chs,'TBSLBRM')
            Nref = floor(3*chs.TBSLBRM/4);
            Ncb = min(N,Nref);
        else % No limit on buffer size
            Ncb = N; 
        end

        info(icw).C = C;                  % Number of code block segements
        info(icw).NLayers = chs.NLayers;  % Number of layers for rate matching
        info(icw).Qm = Qm;                % Modulation order
        info(icw).RV = chs.RV(icw);       % Redundancy version
        
        info(icw).NBG = nbg;              % Base graph number
        info(icw).L = L;                  % Transport block CRC bit length
        info(icw).Kd = Kd;                % Number of bits in each code block (may include CB-CRC bits but excluding Filler bits)
        info(icw).K = K;                  % Number of bits in each code block (including Filler bits)
        info(icw).F = K-Kd;               % Number of fillers bits per code block
        info(icw).N = N;                  % Length of LPDC encoded circular buffer 
        info(icw).E = E;                  % Rate matching output size for each code blocks
        info(icw).Ncb = Ncb;              % Code block soft buffer size
        info(icw).ZC = Zc;                % Selected lifting size
        
    end

end