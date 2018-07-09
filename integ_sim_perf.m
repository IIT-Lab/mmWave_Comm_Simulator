clear; close all; clc;

nUser = [2 4];
UserExp_TRD = load('./SIM_PERF/USER_PER_EXP_TRD.mat', 'PER_LCMV', 'PER_CBF');
UserExp_HEU = load('./SIM_PERF/USER_PER_EXP_HEU.mat', 'PER_HEU');

nAnt = [16 32];
AntExp_TRD = load('./SIM_PERF/ANTENNA_PER_EXP_TRD.mat', 'PER_LCMV', 'PER_CBF');
AntExp_HEU = load('./SIM_PERF/ANTENNA_PER_EXP_HEU.mat', 'PER_HEU');


UserExp_data = [UserExp_TRD.PER_CBF; UserExp_TRD.PER_LCMV; UserExp_HEU.PER_HEU]';
AntExp_data = [AntExp_TRD.PER_CBF; AntExp_TRD.PER_LCMV; AntExp_HEU.PER_HEU]';

figure(1)
bar(nUser, UserExp_data(1:2, :))
legend('CBF', 'LCMV', 'HEU');
title('#User vs PER, 64 antenna')

figure(2)
bar(nAnt, AntExp_data(1:2, :))
legend('CBF', 'LCMV', 'HEU');
title('#Antenna vs PER, 4 users')