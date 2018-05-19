%h5gSSBurst Synchronization Signal burst (SS burst)
%   [WAVEFORM,GRID,INFO] = h5gSSBurst(GNB,BURST) creates an SS burst
%   waveform WAVEFORM and subcarrier grid GRID given general link settings
%   structure GNB and burst configuration structure BURST. Information on
%   the burst structure and content is returned in the structure INFO. The
%   burst is generated according to TS 38.213 Section 4.1 and the SS/PBCH
%   blocks within the burst are generated according to TS 38.211 Section
%   7.4.3.
%
%   Example:
%   % Generate an SS burst with 30kHz subcarrier spacing (Case B), and a
%   % offset of 140 subcarriers from the reference point, defined here as
%   % first subcarrier of the first PRB of both the data numerology (given 
%   % by gnb.SubcarrierSpacing) and the SS/PBCH block numerology (derived 
%   % from burst.BurstType). The subcarrier offset is w.r.t. the subcarrier
%   % spacing of the SS/PBCH block. The SS burst waveform is generated
%   % using a resource grid defined with an internally calculated NDLRB
%   % value that supports the common reference point betweeen data and
%   % SS/PBCH block numerologies and creates a waveform with the same
%   % sampling rate as that of the data numerology (i.e. that given by 
%   % h5gOFDMInfo(gnb)).
%
%   gnb.NDLRB = 216;
%   gnb.SubcarrierSpacing = 15;     % 15kHz data SCS
%   gnb.NCellID = 42;
%
%   burst.BurstType = 'CaseB';      % 30kHz SS/PBCH block SCS
%   burst.DisplayBurst = true;
%   burst.SubcarrierOffset = 140;
%   burst.SSBTransmitted = [1 1 1 0 1 1 1 1];
%
%   [ssbWaveform,ssbGrid,ssbInfo] = h5gSSBurst(gnb,burst);
%
%   See also h5gPSS, h5gPSSIndices, h5gSSS, h5gSSSIndices, h5gPBCH,
%   h5gPBCHIndices, h5gPBCHDMRS, h5gPBCHDMRSIndices.

%   Copyright 2018 The MathWorks, Inc.

function [waveform,grid,info] = h5gSSBurst(gnb,burst)
    
    if (~isfield(burst,'BurstType'))
        burst.BurstType = 'CaseA';
    end
    
    if (~isfield(burst,'SSBTransmitted') || isempty(burst.SSBTransmitted))
        % Default SSB-transmitted bitmap length from choices in TS 38.331
        % ServingCellConfigCommon IE ssb-PositionsInBurst
        if (any(strcmpi(burst.BurstType,{'CaseA' 'CaseB' 'CaseC'})))
            % Short bitmap, for sub 3 GHz
            burst.SSBTransmitted = ones(1,4);
        else
            % Long bitmap, for above 6GHz 
            burst.SSBTransmitted = ones(1,64);
        end
    end
    
    % Get starting symbols of SS blocks in the half frame according to TS
    % 38.213 Section 4.1. The number of SS blocks is denoted L
    ssbStartSymbols = getStartSymbols(burst);
    L = length(ssbStartSymbols);
    
    % Get the half frame number of the SS burst
    if (L==4)
        % Transmission in half frame 1 of a frame (the 2nd half frame) is
        % only applicable for L==4 and 5ms SSB periodicity
        if (~isfield(burst,'NHalfFrame'))
            burst.NHalfFrame = 0;
        end
        n_hf = mod(burst.NHalfFrame,2);
    else
        n_hf = 0;
    end

    % Get OFDM information for data numerology
    dataOFDMInfo = h5gOFDMInfo(gnb);
    
    % Get subcarrier spacing for SS burst
    gnb.SubcarrierSpacing = getSSBSubcarrierSpacing(burst);
    
    % Calculate number of subcarriers in SS burst numerology that span the
    % subcarriers in the data numerology, and calculate number of resource
    % blocks required to span that number of subcarriers
    ssbNSubcarriers = dataOFDMInfo.NSubcarriers * dataOFDMInfo.SubcarrierSpacing / gnb.SubcarrierSpacing;
    ssbNDLRB = ceil(ssbNSubcarriers / 12);
    
    if (ssbNDLRB < 20)
        error('lte:error','Number of resource blocks (%d) in SS/PBCH block subcarrier spacing (%dkHz) must be at least 20.',ssbNDLRB,gnb.SubcarrierSpacing);
    end
    
    % Create resource grid for SS burst
    gnb.NDLRB = ssbNDLRB;
    ssbOFDMInfo = h5gOFDMInfo(gnb);
    grid = zeros([ssbOFDMInfo.NSubcarriers ssbOFDMInfo.SymbolsPerSubframe*5 1]);
    displaygrid = grid;
    
    % Calculate offset to the subcarrier in SS burst grid which represents
    % data subcarrier k=0, which will be used as the reference point
    % between numerologies
    ref_point_offset = (ssbOFDMInfo.NSubcarriers - ssbNSubcarriers) / 2;
    
    % Default the subcarrier offset to place SS burst in the centre of the
    % band used by the data
    if (~isfield(burst,'SubcarrierOffset'))
        burst.SubcarrierOffset = ssbOFDMInfo.NSubcarriers/2 - 120;
    end
    
    % Limit the subcarrier offset if the value places the SS block beyond
    % the upper edge of the configured bandwidth
    if ((ref_point_offset + burst.SubcarrierOffset + 240) > size(grid,1))
        burst.SubcarrierOffset = size(grid,1) - 240 - ref_point_offset;
    end

    % For each transmitted SSB
    for i_SSB = find(burst.SSBTransmitted) - 1

        ssbGrid = zeros([240 4 1]);

        pssInd = h5gPSSIndices();
        pss = h5gPSS(gnb);
        ssbGrid(pssInd) = pss;

        sssInd = h5gSSSIndices();
        sss = h5gSSS(gnb);
        ssbGrid(sssInd) = sss;
        
        [pbchInd,pbchIndInfo] = h5gPBCHIndices(gnb);
        pn = comm.PNSequence;
        pn.SamplesPerFrame = pbchIndInfo.G;
        cw = pn();
        pbch = h5gPBCH(gnb,i_SSB,cw);
        ssbGrid(pbchInd) = pbch;

        pbchDmrsInd = h5gPBCHDMRSIndices(gnb);
        pbchDmrs = h5gPBCHDMRS(gnb,4*i_SSB + n_hf);
        ssbGrid(pbchDmrsInd) = pbchDmrs;

        grid(ref_point_offset + burst.SubcarrierOffset + (1:240),ssbStartSymbols(i_SSB+1)+(1:size(ssbGrid,2)),:) = ssbGrid;

        if (isfield(burst,'DisplayBurst') && burst.DisplayBurst)
        
            displayburst = zeros([240 4 1]);
            displayburst(pssInd) = 1;
            displayburst(sssInd) = 2;
            displayburst(pbchInd) = 3;
            displayburst(pbchDmrsInd) = 4;
            displaygrid(ref_point_offset + burst.SubcarrierOffset + (1:240),ssbStartSymbols(i_SSB+1)+(1:size(displayburst,2)),:) = displayburst;
            
        end
        
    end

    % Modulate grid
    waveform = h5gOFDMModulate(gnb,grid);
    
    % Display burst if requested
    if (isfield(burst,'DisplayBurst') && burst.DisplayBurst)
        figure;
        imagesc(displaygrid);
        axis xy;
        ylabel('Subcarriers');
        xlabel('OFDM symbols');
        title(sprintf('SS burst, SCS=%dkHz, NDLRB=%d',gnb.SubcarrierSpacing,gnb.NDLRB));
        hold on;
        for i = 1:4
            patch([-2 -3 -3 -2],[-2 -2 -3 -3],i);
        end
        legend('PSS','SSS','PBCH','PBCH DM-RS');
        drawnow;
    end
    
    % Create information output
    info = struct();
    info.SamplingRate = ssbOFDMInfo.SamplingRate;
    info.Nfft = ssbOFDMInfo.Nfft;
    info.NSubcarriers = gnb.NDLRB*12;
    info.SubcarrierSpacing = gnb.SubcarrierSpacing;
    info.OccupiedSubcarriers = ref_point_offset + burst.SubcarrierOffset + (0:239).';
    info.Symbols = (0:3).' + ssbStartSymbols;
    info.Symbols = info.Symbols(:).';
    info.SSBTransmitted = burst.SSBTransmitted;
    info.SymbolsTransmitted = (0:3).' + ssbStartSymbols(logical(burst.SSBTransmitted));
    info.SymbolsTransmitted = info.SymbolsTransmitted(:).';    
    
end

function scs = getSSBSubcarrierSpacing(burst)

    if (strcmpi(burst.BurstType,'CaseA'))
        scs = 15;
    elseif (any(strcmpi(burst.BurstType,{'CaseB','CaseC'})))
        scs = 30;
    elseif (strcmpi(burst.BurstType,'CaseD'))
        scs = 120;
    elseif (strcmpi(burst.BurstType,'CaseE'))
        scs = 240;
    end

end

function ssbStartSymbols = getStartSymbols(burst)

    % 'alln' gives the overall set of SS block indices 'n' described in 
    % TS 38.213 Section 4.1, from which a subset is used for each Case A-E    
    alln = [0; 1; 2; 3; 5; 6; 7; 8; 10; 11; 12; 13; 15; 16; 17; 18];
    
    L = length(burst.SSBTransmitted);
    
    switch (burst.BurstType)
        case 'CaseA'
            m = 14;
            i = [2 8];
            if (L==4)
                nn = 2;
            elseif (L==8)
                nn = 4;
            else
                error('For Case A, the SSBTransmitted bitmap must be of length 4 or 8.');
            end
        case 'CaseB'
            m = 28;
            i = [4 8 16 20];
            if (L==4)
                nn = 1;
            elseif (L==8)
                nn = 2;
            else
                error('For Case B, the SSBTransmitted bitmap must be of length 4 or 8.');
            end
        case 'CaseC'
            m = 14;
            i = [2 8];
            if (L==4)
                nn = 2;
            elseif (L==8)
                nn = 4;
            else
                error('For Case C, the SSBTransmitted bitmap must be of length 4 or 8.');
            end
        case 'CaseD'
            m = 28;
            i = [4 8 16 20];
            if (L==64)
                nn = 16;
            else
                error ('For Case D, the SSBTransmitted bitmap must be of length 64.');
            end
        case 'CaseE'
            m = 56;
            i = [8 12 16 20 32 36 40 44];
            if (L==64)
                nn = 8;
            else
                error ('For Case E, the SSBTransmitted bitmap must be of length 64.');
            end
    end
    
    n = alln(1:nn);
    ssbStartSymbols = (i + m*n).';
    ssbStartSymbols = ssbStartSymbols(:).';
    
end
