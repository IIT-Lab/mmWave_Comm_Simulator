clear; close all;

load('paul-test.mat');
W = W';
%% Parameters privite to this function
SNR = 10;

%% All parameters needed from upper layer
% PHY
lengthPSDU = 1000;
MCS = 1;

% PHASED
nTx_row = 8;
nTx_col = 8;
nRx = 1;
[nAntennas, nUsers] = size(W);

% Simulation scenarios
angleToRx = zeros(nUsers, 1); %randi([-30, 30], nUsers, 1); % Shall be [Azimuth; Elevation];
distance = 1 + rand([nUsers, 1]);

%% Configure system objects
tx_phy = s_phy_tx( ...
    'MCS', MCS, ...
    'PSDULength', lengthPSDU);
tx_pha = s_phased_tx( ...
    'numTxElements_row', nTx_row, ...
    'numTxElements_col', nTx_col, ...
    'txGain',            20, ...
    'visualization', false);
channel = s_phased_channel( ...
    'numInputElements_row',     nTx_row, ...
    'numInputElements_col',     nTx_col, ...
    'numOutputElements_row',    nRx, ...
    'numOutputElements_col',    nRx, ...
    'SNR',                      SNR, ...
    'applyPathLoss',                  true);
rx_pha = s_phased_rx( ...
    'numRxElements', nRx, ...
    'rxGain',        20);
rx_phy = s_phy_rx();

%% Simulations
psdu = cell(nUsers, 1);
txWaveforms = cell(nUsers, 1);
response = cell(nUsers, 1);
go_through = false(nUsers, 1);

for user = 1 : nUsers
    psdu{user} = randi([0 1], lengthPSDU * 8, 1);
    [txSymbols, cfgDMG] = tx_phy(psdu{user});
    [txWaveforms{user}, response{user}] = tx_pha(txSymbols, angleToRx(user), W(:, user));
end

for outer_iter = 1 : nUsers
    combined_tx_waveforms = txWaveforms{outer_iter};
    for inner_iter = 1 : nUsers
        if inner_iter ~= outer_iter
            combined_tx_waveforms = combined_tx_waveforms + txWaveforms{inner_iter} * response{outer_iter};
        end
    end

    txWaveforms_afterChannel = channel(combined_tx_waveforms, distance(user));
    rxSymbols = rx_pha(txWaveforms_afterChannel, angleToRx(outer_iter));
    [psdu_rx, rxflag] = rx_phy(rxSymbols, cfgDMG);

    if ~isempty(psdu_rx)
        go_through(outer_iter) = ~(any(biterr(psdu{outer_iter}, psdu_rx)) && rxflag);
    end
end

%% BER sim
% totPkt = 1;
% for user = 1 : nUsers
%     psdu = randi([0 1], lengthPSDU * 8, 1);
%     [txSymbols, cfgDMG] = tx_phy(psdu);
%     [txWaveforms, response] = tx_pha(txSymbols, angleToRx(user), []);
%     txWaveforms_afterChannel = channel_pha(txWaveforms, distance(user));
%     rxSymbols = rx_pha(txWaveforms_afterChannel, angleToRx(user));
%     [psdu_rx, rxflag] = rx_phy(rxSymbols, cfgDMG);
%     if ~isempty(psdu_rx)
%         numErr = any(biterr(psdu,psdu_rx))
%     end
% end