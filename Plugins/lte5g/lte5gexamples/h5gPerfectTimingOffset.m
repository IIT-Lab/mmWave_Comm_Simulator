%h5gPerfectTimingOffset perfect timing offset

% Copyright 2016-2017 The MathWorks, Inc.

function [offset,mag] = h5gPerfectTimingOffset(pathGains,channelInfo,SR)
    
    % Get number of paths 'L', number of transmit antennas 'P' and number
    % of receive antennas 'R' in the path gains array
    [~,L,P,R] = size(pathGains);
    
    % Calculate path sample delays from path delays and sampling rate
    pathSampleDelays = channelInfo.PathDelays * SR;
    
    % Establish integer path sample indices 'n'
    n = floor(pathSampleDelays) + 1;
    
    % Establish fractional path delays 'rho'
    rho = pathSampleDelays + 1 - n;
    
    % Create channel impulse response array 'h' for each impulse response
    % sample, transmit antenna and receive antenna
    h = zeros(max(n)+1,P,R);
    
    % Average the path gains array across all time elements (1st dimension)
    pathGains = permute(mean(pathGains,1),[2 3 4 1]);
    
    % For each path, add its contribution to the channel impulse response
    % across all transmit and receive antennas. This involves positioning
    % the path according to its integer path sample index n(l) and
    % performing a first-order interpolation between two sample indices
    % according to the fractional path delay rho(l)
    for l = 1:L
        
        k = n(l) + [0;1];
        h(k,:,:) = h(k,:,:) + [(1-rho(l)); rho(l)] .* pathGains(l,:,:);
        
    end
    
    % Combine the transmit antennas in the channel impulse response array,
    % leaving a matrix of impulse response samples versus receive antennas
    h = permute(sum(h,2),[1 3 2]);
    
    % Take the magnitude of the impulse response matrix
    mag = abs(h);
    
    % Return the minimum timing offset across the receive antennas that
    % have a correlation peak at least 50% of the magnitude of the
    % strongest peak, including accounting for the channel filter
    % implementation delay
    offset = zeros(1,R);
    maxmag = zeros(1,R);
    for r = 1:R
        maxmag(r) = max(mag(:,r));
        offset(r) = find(mag(:,r)==maxmag(r),1) - 1;
    end
    offset = min(offset(maxmag>=0.5*max(maxmag)));
    offset = offset + channelInfo.ChannelFilterDelay;
    
end
