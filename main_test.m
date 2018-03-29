clear; close all;

lengthPSDU = 1000;
nTx_row = 8;
nTx_col = 8;

nRx = 1;

tic
for SNR = 5
    for MCS = 1
        tx_phy = s_phy_tx( ...
            'MCS', MCS, ...
            'PSDULength', lengthPSDU);
        tx_pha = s_phased_tx( ...
            'numTxElements_row', nTx_row, ...
            'numTxElements_col', nTx_col, ...
            'visualization', false);
        channel_pha = s_phased_channel( ...
            'numInputElements_row',     nTx_row, ...
            'numInputElements_col',     nTx_col, ...
            'numOutputElements_row',    1, ...
            'numOutputElements_col',    1, ...
            'SNR',                      SNR);
        rx_pha = s_phased_rx( ...
            'numRxElements', nRx);
        
        rx_phy = s_phy_rx();
        
        totPkt = 100;
        numErr = 0;
        
        angleToRx = randi([-30, 30]);
        
        for i = 1 : totPkt
            psdu = randi([0 1], lengthPSDU*8, 1);
            [txSymbols, cfgDMG] = tx_phy(psdu);
            txWaveforms = tx_pha(txSymbols, angleToRx, []);
            txWaveforms_afterChannel = channel_pha(txWaveforms);
            rxSymbols = rx_pha(txWaveforms_afterChannel, angleToRx);
            [psdu_rx, rxflag] = rx_phy(rxSymbols, [], cfgDMG);
            if ~isempty(psdu_rx)
                numErr = any(biterr(psdu,psdu_rx)) + numErr;
            end
        end
        fprintf('MCS = %d, SNR = %.2fdB, PER = %.2f%%\n', MCS, SNR, numErr * 100 / totPkt);
    end
end
toc