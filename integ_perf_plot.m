clc; clear all; close all;

load('./SIM_PERF/Final/USER_PER_HEU_128_ANT_1_MCS.mat');
figure
bar(PER_HEU)
legend('2 Users', '4 Users', '8 Users')
title(['PER performance of Heuristics, ' num2str(nAntenna) ' antennas'])
xticklabels(cdlProfile)

load('./SIM_PERF/Final/USER_PER_TRD_128_ANT_1_MCS.mat');
figure
bar(PER_LCMV)
legend('2 Users', '4 Users', '8 Users')
title(['PER performance of LCMV, ' num2str(nAntenna) ' antennas'])
xticklabels(cdlProfile)

load('./SIM_PERF/Final/ANT_PER_HEU_4_USERS_1_MCS.mat');
figure
bar(PER_HEU)
legend('16 antennas', '32 antennas', '64 antennas', '128 antennas', '144 antennas')
title(['PER performance of Heuristics, ' num2str(nUsers) ' users'])
xticklabels(cdlProfile)

load('./SIM_PERF/Final/ANT_PER_TRD_4_USERS_1_MCS.mat');
figure
bar(PER_LCMV)
legend('16 antennas', '32 antennas', '64 antennas', '128 antennas', '144 antennas')
title(['PER performance of LCMV, ' num2str(nUsers) ' users'])
xticklabels(cdlProfile)