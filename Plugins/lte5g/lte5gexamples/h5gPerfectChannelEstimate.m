%h5gPerfectChannelEstimate perfect channel estimation

% Copyright 2016-2018 The MathWorks, Inc.

function H = h5gPerfectChannelEstimate(enb,channel,pathGains,toffset)
    
    % Establish number of subframes 'nsf' for which to create estimate
    if isfield(enb,'TotSubframes')
        nsf = enb.TotSubframes;
    else
        nsf = 1;
    end
    
    % Get number of transmit antennas 'P' and number of receive antennas
    % 'R' in the path gains array
    [~,~,P,R] = size(pathGains); 
    
    % Create waveform containing an impulse at the start of each OFDM
    % symbol, for a single transmit antenna
    ofdmInfo = h5gOFDMInfo(enb);
    txGrid = ones(ofdmInfo.Nfft,ofdmInfo.SymbolsPerSlot*nsf,1);
    txWave = h5gOFDMModulate(enb,txGrid,0);
    txWave = [txWave; txWave(1:toffset,:)];
    
    % Establish which input waveform sample times correspond to which
    % channel coefficient sample times. 'idx' is a vector indicating the
    % 1st dimension index of 'pathGains' for each sample of the input
    % waveform 'txWave'
    if (isa(channel,'nr5gTDLChannel'))
        
        pathGains = pathGains(1:size(txWave,1),:,:,:);
        idx = (1:size(pathGains,1)).';
        
    elseif (isa(channel,'nr5gCDLChannel'))
        
        idx = getSampleIndices(channel,pathGains,size(txWave));
        
    end

    % Get filter coefficients for each path filter
    pathFilters = getPathFilters(channel);
    
    % Filter waveform by each path filter giving a matrix of impulse
    % response samples versus path
    delayedWaves = conv2(txWave(1:end-size(pathFilters,1)+1),pathFilters);
    
    % Create the received waveform for each transmit and receive antenna by
    % duplicating the filtered impulse responses for each transmit and
    % receive antenna, applying the path gains and combining the paths
    delayedWaves = repmat(delayedWaves,[1 1 P]);   % expand transmit antennas. size: T-by-L-by-P
    delayedWaves = repmat(delayedWaves,[1 1 1 R]); % expand receive antennas. size: T-by-L-by-P-by-R
    rxWave = delayedWaves .* pathGains(idx,:,:,:); % apply cluster gains. size: T-by-L-by-P-by-R
    rxWave = sum(permute(rxWave,[1 3 4 2]),4);     % combine path. size: T-by-P-by-R

    % Perform synchronization
    rxWave = rxWave(1+toffset:end,:,:);

    % For each receive antenna, OFDM demodulate the received waveform
    % across all transmit antennas to form the overall channel estimate
    % array
    H = zeros(ofdmInfo.NSubcarriers,size(txGrid,2),R,P);
    for idx = 1:R
        
        H(:,:,idx,:)= h5gOFDMDemodulate(enb,rxWave(:,:,idx));
        
    end
    
end

function idx = getSampleIndices(channel,pathGains,insize)

    if (channel.MaximumDopplerShift && channel.SampleDensity)
        if (isfinite(channel.SampleDensity))
            pathSampleRate = channel.MaximumDopplerShift * 2 * channel.SampleDensity;
        else
            pathSampleRate = channel.SampleRate;
        end
        Ht = (0:(size(pathGains,1)-1)).'/pathSampleRate;
    else
        Ht = 0;
    end
    t = (0:(insize(1)-1)).'/channel.SampleRate;
    if (numel(Ht)>1)
        Ht = Ht + mean(diff(Ht))/2;
    end
    T_H = numel(Ht);
    if (~isempty(t))
        if (t(1) < Ht(1))
            Ht = [t(1); Ht];
        end
    end
    idx = sum(Ht.'<=t,2);
    idx(idx>T_H) = T_H;
    
end
    