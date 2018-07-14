clc; clear all; close all;

load('../SIM_PERF/USER_PER_EXP_TRD_64_MCS_1.mat');

figure
bar(PER_LCMV)
xticklabels(cdlProfile)
legend({'2 Users', '4 Users', '8 Users'})
title(['PER of LCMV, ' num2str(nAntenna) ' antenna'])

figure
bar(PER_CBF)
xticklabels(cdlProfile)
legend({'2 Users', '4 Users', '8 Users'})
title(['PER of CBF, ' num2str(nAntenna) ' antenna'])

load('../SIM_PERF/USER_PER_EXP_TRD_144_MCS_1.mat');

figure
bar(PER_LCMV)
xticklabels(cdlProfile)
legend({'2 Users', '3 Users', '4 Users', '6 Users'})
title(['PER of LCMV, ' num2str(nAntenna) ' antenna'])

figure
bar(PER_CBF)
xticklabels(cdlProfile)
legend({'2 Users', '3 Users', '4 Users', '6 Users'})
title(['PER of CBF, ' num2str(nAntenna) ' antenna'])

