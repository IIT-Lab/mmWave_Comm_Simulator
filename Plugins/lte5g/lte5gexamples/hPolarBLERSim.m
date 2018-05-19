function [numFerr,nFrame] = hPolarBLERSim(K,E,crcLen,F,L,EbNo, ...
                                maxNumFrames,maxNumErrs,minFrames)
%Polar Coding BLER Simulator
%
%   Polar Code Construction and Encoding and Polar Decoding, using the
%   CRC-Aided Successive Cancellation LIST Decoding algorithm, for AWGN
%   channel.
%
%   See also NewRadio5GPolarCodingExample, hPolarBLERSimDriver.

%   Copyright 2017-2018 The MathWorks, Inc.

%#codegen

persistent crcGen polarEnc bpskMod chan bpskDemod polarDec;

% Downlink scenario
nPC = 0;            %#ok  number of parity-check bits
nMax = 9;           %#ok  maximum n for downlink
iIL = true;         % interleave input
iBIL = false;       % interleave coded bits

N = length(F);
R = K/E;

if isempty(crcGen)
    if crcLen==6
        crcPoly = [1 1 0 0 0 0 1];
    elseif crcLen==11
        crcPoly = [1 1 1 0 0 0 1 0 0 0 0 1];
    else % 24C
        crcPoly = [1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1];
    end
    crcGen = comm.CRCGenerator('Polynomial', crcPoly);
end

if isempty(polarEnc)
    polarEnc = h5gPolarEncoder(N,K,F,'InterleaveInput',iIL); 
end

if isempty(bpskMod)
    bpskMod = comm.BPSKModulator;
end

if isempty(chan)
    % AWGN Channel
    chan = comm.AWGNChannel('NoiseMethod','Variance', ...
        'RandomStream','mt19937ar with seed');
end

if isempty(bpskDemod)
    % BPSK Demodulator
    bpskDemod = comm.BPSKDemodulator( ...
        'DecisionMethod', 'Approximate log-likelihood ratio');
end

if isempty(polarDec)
    polarDec = h5gPolarDecoder(N,K,F,L,crcLen,'DeinterleaveOutput',iIL); 
end

bps = 1;                          % bits per symbol, 1 for BPSK, 2 for QPSK
EsNo = EbNo + 10*log10(bps);       
snrdB = EsNo + 10*log10(R);       % in dB
noiseVar = 1./(10.^(snrdB/10)); 

chan.Variance = noiseVar;
bpskDemod.Variance = noiseVar;

numFerr = 0; nFrame = 0;
while (nFrame < maxNumFrames) && (numFerr < maxNumErrs) ... 
        || (nFrame < minFrames)

    % Generate a random message
    msg = randi([0 1],K-crcLen,1);

    % CRC encode
    msgcrc = crcGen(msg);

    % Polar encode
    encOut = polarEnc(msgcrc);

    % Rate match
    modIn = h5gRateMatchPolar(encOut,K,E,iBIL);
    
    % Modulate
    chIn = bpskMod(modIn);

    % Add WGN
    rSig = chan(chIn);

    % Soft demodulate
    rxLLR = bpskDemod(rSig);

    % Rate recover
    decIn = h5gRateRecoverPolar(rxLLR,K,N,iBIL);
    
    % Polar decode
    decBits = polarDec(decIn);

    % Compare
    numFerr = numFerr + any(decBits(1:K-crcLen)~=msg);
    nFrame = nFrame + 1;
end

end
