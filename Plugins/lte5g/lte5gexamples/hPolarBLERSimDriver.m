%% Polar Coding BLER Simulation Driver
%   Uses both MATLAB Coder and Parallel Computing Toolbox for acceleration.
%
%   See also NewRadio5GPolarCodingExample, hPolarBLERSim.

%   Copyright 2017-2018 The MathWorks, Inc.

%s = rng(0);       % Seed the RNG for repeatability
%myPool = parpool('AttachedFiles', {'hPolarBLERSim_mex.mexw64'});

%% 
% Code parameters
K = 54;        
E = 124;     

% EbNo ranges for rates (EbNo==EsNo for BPSK)
%     EbNoVec = -3:0.5:5.5;   % for r~=1/16(56/864)
%     EbNoVec = -5:0.5:3.5;   % for r~=1/4 (43/180)
%     EbNoVec = -4:0.5:4.5;   % for r~=2/5 (54/124)
%     EbNoVec = -3:0.5:5.5;   % for r~=2/3 (164/240)
%     EbNoVec = 0:0.5:7.5;    % for r~=8/9 (164/184)
EbNoVec = -4:0.5:4.5;  

% Parameter sets based on different rates
% E = 864; K = 56;  % N = 512, rep, K/E~=1/16 
% E = 180; K = 43;  % N = 256, punc <, K/E~=1/4
% E = 124; K = 54;  % N = 128, punc >=,  K/E~=2/5 
% E = 240; K = 164; % N = 256, short, K/E~=2/3 
% E = 184; K = 164; % N = 256, short, K/E~=8/9 
     
%%
% Downlink scenario (K > 24, includes crc)
crcLen = 24;      % Number of CRC bits for DL
nPC = 0;          % Number of parity-check bits, assumed=0
nMax = 9;         % Maximum n value for 2^n
iIL = true;       % Interleave input
iBIL = false;     % Interleave coded bits
L = 8;            % List length, a power of two, [2 4 8]

% Uplink scenario (K > 30, includes crc)
% crcLen = 11;    % Number of CRC bits for UL: 11 for K>30, 6 for 18<=K<=25
% nPC = 0;        % for K > 30, 3 for 18<=K<=25
% nMax = 10;
% iIL = false; 
% iBIL = true;    % not implemented

F = h5gPolarConstruct(K,E,nMax); % assumes nPC = 0 only
N = length(F);

%%
% Simulation parameters
genCode = true;     % Enable codegen

maxNumErrsVec = 1000*ones(1,length(EbNoVec)); % errors collected at each EbNo pt
% maxNumFrames = 100e3;  % maximum frames simulated for each EbNo pt
% minFrames = 1000;      % minimum frames run, for low SNR range
maxNumFrames = 100;      % maximum frames simulated for each EbNo pt
minFrames = 100;         % minimum frames run, for low SNR range

%% Generate code
if genCode
    disp('Generating code')

    % Code generation parameters
    cfg = coder.config;
    cfg.IntegrityChecks = false;
    cfg.ResponsivenessChecks = false;
    cfg.EnableJIT = true;
    cfg.EchoExpressions = false;
    cfg.CompileTimeRecursionLimit = 4096;
    cfg.ExtrinsicCalls = false;

    codegen('hPolarBLERSim', '-args', ...
        {coder.Constant(K),coder.Constant(E),coder.Constant(crcLen), ...
         coder.Constant(F),coder.Constant(L),EbNoVec(1),maxNumFrames, ...
         maxNumErrsVec(1),minFrames}, '-config',cfg);
end

%% Processing loop
tic;
bler = zeros(1,length(EbNoVec)); nferr = bler; nFrame = bler;
parfor idx = 1:length(EbNoVec)
    
    EbNo = EbNoVec(idx);
    disp(['Processing K: ' num2str(K)  ', E: ' num2str(E) ', EbNo: ' ...
          num2str(EbNo) ' dB'])
    
    % Run simulation
    if genCode
        [numFerr, nF] = hPolarBLERSim_mex(K,E,crcLen,F,L,EbNo, ...
            maxNumFrames,maxNumErrsVec(idx),minFrames);
    else        
        [numFerr, nF] = hPolarBLERSim(K,E,crcLen,F,L,EbNo, ... 
            maxNumFrames,maxNumErrsVec(idx),minFrames);   %#ok
    end
    
    bler(1,idx) = numFerr/nF;
    nferr(1,idx) = numFerr;
    nFrame(1,idx) = nF;
end
toc;