% The purpose of this script is to simulate 802.11ad performance of MCS 1 - 12

addpath('./Plugins/lte5g');
% diary('DELAY_SPREAD.log');
lengthPSDU = 1024;

nTx_row = 1;
nTx_col = 1;
nRx = 1;

tic
for profile = {'CDL-A', 'CDL-B', 'CDL-C', 'CDL-D', 'CDL-E'}
    for SNR = 0 : 20
        for MCS = 1 : 12
            tx_phy = s_phy_tx( ...
                'MCS', MCS, ...
                'PSDULength', lengthPSDU);
            channel_pha = s_phased_channel( ...
                'numInputElements_row',     1, ...
                'numInputElements_col',     1, ...
                'numOutputElements_row',    1, ...
                'numOutputElements_col',    1, ...
                'SNR',                      SNR, ...
                'applyPathLoss',            false, ...
                'profile',                  cell2mat(profile), ...
                'seed',                     randi([1000, 10000]));
            rx_phy = s_phy_rx();
            
            totPkt = 1;
            numError = 0;
            
            for i = 1 : totPkt
                psdu = randi([0 1], lengthPSDU * 8, 1);
                [txSymbols, cfgDMG] = tx_phy(psdu);
                txWaveforms_afterChannel = channel_pha(txSymbols, 0, 0);
                [psdu_rx, rxflag] = rx_phy(txWaveforms_afterChannel, cfgDMG);
                
                if ~isempty(psdu_rx)
                    bitErrorFlag = any(biterr(psdu, psdu_rx));
                else
                    bitErrorFlag = true;
                end
                % if we receive nothing, or rxFlag is false, 
                % or we have a bit error, we claim packet error
                if isempty(psdu_rx) || ~rxflag || bitErrorFlag 
                    numError = numError + 1;
                end
            end
            fprintf('%s\t%d\t%d\t%.2f\n', cell2mat(profile), MCS, SNR, numError * 100 / totPkt);
        end
    end
end
toc
% diary off