clear; close all;

load('paul-test.mat');
W = W';
%% All parameters needed
lengthPSDU = 1000;
nTx_row = 8;
nTx_col = 8;
nRx = 1;
SNR = 10;
MCS = 1;
[nAntennas, nUsers] = size(W);
angleToRx = randi([-30, 30], nUsers, 1);
distance = 1 + rand([nUsers, 1]);

%% Configure system objects
tx_phy = s_phy_tx( ...
    'MCS', MCS, ...
    'PSDULength', lengthPSDU);
tx_pha = s_phased_tx( ...
    'numTxElements_row', nTx_row, ...
    'numTxElements_col', nTx_col, ...
    'txGain',        1, ...
    'visualization', false);
channel_pha = s_phased_channel( ...
    'numInputElements_row',     nTx_row, ...
    'numInputElements_col',     nTx_col, ...
    'numOutputElements_row',    nRx, ...
    'numOutputElements_col',    nRx, ...
    'SNR',                      SNR, ...
    'applyPathLoss',            false);
rx_pha = s_phased_rx( ...
    'numRxElements', nRx, ...
    'rxGain',        1);
rx_phy = s_phy_rx();

%% Simulations
psdu = cell(nUsers, 1);
txWaveforms = cell(nUsers, 1);
response = cell(nUsers, 1);
go_through = false(nUsers, 1);


for user = 1 : nUsers
    psdu{user} = randi([0 1], lengthPSDU * 8, 1);
    [txSymbols, cfgDMG] = tx_phy(psdu{user});
    [txWaveforms{user}, response{user}] = tx_pha(txSymbols, angleToRx(user), []);%W(:, user));
end

for user_outer = 1 : nUsers
    combined_tx_waveforms = txWaveforms{user_outer};
    for user_inner = 1 : nUsers
        combined_tx_waveforms = combined_tx_waveforms + txWaveforms{user_inner} * 1;%response{user_outer};
    end

    txWaveforms_afterChannel = channel_pha(combined_tx_waveforms, distance(user));
    rxSymbols = rx_pha(txWaveforms_afterChannel, angleToRx(user_outer));
    [psdu_rx, rxflag] = rx_phy(rxSymbols, cfgDMG);

    if ~isempty(psdu_rx)
        go_through(user_outer) = ~(any(biterr(psdu{user_outer}, psdu_rx)) && rxflag);
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